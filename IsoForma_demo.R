
# Tested on a Windows 10 machine
# R version 4.3.1

# Install packages
install.packages("devtools") 
devtools::install_github("EMSL-Computing/pspecterlib")
devtools::install_github("EMSL-Computing/isoforma-lib")
install.packages("xlsx")
install.packages("doParallel")

# Load the packages
library(pspecterlib)
library(isoforma)
library(xlsx)
library(doParallel)

# Set path to files and read list of targets and parameters:
setwd("E:/IsoForma/IsoForma-paper") # set as base path, files Targets.xlsx and Parameters.xlsx must be here
mzMLFolder = "RawData"
targets = read.xlsx("Targets.xlsx", 1) # Targets file must have the following columns: "ProteoformGroup", "RTStart", "RTEnd", "UnmodifiedSequence", "Modifications".
parameters = read.xlsx("Parameters.xlsx", 1) 

# Set number of cores to run multiple files in parallel:
nCores = 3

#' Run the main IsoForma function
IsoFormaMain <- function(mzMLFullPath, outPath, 
                         UnmodSequence,
                         Modification,
                         RTStart,
                         RTEnd,
                         parameters) {
  # Parse parameters:
  MaxPrecursorCharge = parameters[parameters$Parameter == "MaxPrecursorCharge", "Setting"] %>% as.numeric()
  PPMPrecursor = parameters[parameters$Parameter == "PPMPrecursor", "Setting"] %>% as.numeric()
  MinMS1Matches = parameters[parameters$Parameter == "MinMS1Matches", "Setting"] %>% as.numeric()
  MassWindow = parameters[parameters$Parameter == "MassWindow", "Setting"] %>% as.numeric()
  MaxFragmentCharge =  parameters[parameters$Parameter == "MaxFragmentCharge", "Setting"] %>% as.numeric()
  PPMRound = parameters[parameters$Parameter == "PPMRound", "Setting"] %>% as.numeric()
  PPMThreshold = parameters[parameters$Parameter == "PPMThreshold", "Setting"] %>% as.numeric()
  CorrelationScore = parameters[parameters$Parameter == "CorrelationScore", "Setting"] %>% as.numeric()
  
  # Load mzML file
  ScanMetadata = get_scan_metadata(mzMLFullPath)

  # Create output folder
  if (!dir.exists(outPath)) {dir.create(outPath)}
  
  # Generate the list of PTMs
  MultipleMods = pspecterlib::multiple_modifications(
    Sequence = UnmodSequence,
    Modification = Modification,
    ReturnUnmodified = TRUE
  )

  # 1. Output scan number selection
  ScanNumbers = pull_scan_numbers(Sequence = as.character(MultipleMods[length(MultipleMods)]),
                    ScanMetadata,
                    RTStart,
                    RTEnd,
                    MassWindow = MassWindow,
                    MinMS1Matches = MinMS1Matches,
                    ppmMzTolerance = PPMPrecursor,
                    IsotopeAlgorithm = "isopat",
                    ActivationMethod = "ETD",
                    MinAbundance = 1)
  if (is.null(ScanNumbers)) {return(NULL)}
  #ScanNumbers = c(3924,3932,3936,3949,3951,3961)
  
  selectedScans = ScanMetadata[which(ScanMetadata$`Scan Number` %in% ScanNumbers),
                               c("Scan Number", "MS Level", "Retention Time", "Precursor M/Z", "Precursor Charge", "Precursor Scan", "Activation Method")]
  
  write.csv(selectedScans, file = file.path(outPath, "SelectedScans.csv"), row.names = F, quote = F)
  
  PrecursorsMzs = paste(selectedScans$`Precursor M/Z`, collapse = ' ; ')
  PrecursorsCharges = paste(selectedScans$`Precursor Charge`, collapse = ' ; ')
  ScanNumbers = as.numeric(selectedScans$`Scan Number`)
  
  # Sum MS2 spectra
  ms2 = do.call(rbind, lapply(ScanNumbers, function(scan_number) {
    pspecterlib::get_peak_data(ScanMetadata, scan_number)}))
  ms2 = ms2[with(ms2, order(ms2$`M/Z`)), ] 
  # Fix peak data list classes
  ms2 = list(pspecterlib::make_peak_data(ms2$`M/Z`, ms2$Intensity))
  ms2 = sum_ms2_spectra(PeakDataList = ms2, PPMRound = PPMRound)

  # # Sum ms2 spectra
  # spectra <- sum_ms2_spectra(
  #   ScanMetadata = ScanMetadata2,
  #   ScanNumbers = as.numeric(ScanNumbers),
  #   PPMRound = PPMRound
  # )
  
  # Save summed spectrum as MGF:
  write_mgf_simple(MZ = ms2$`M/Z`, Intensity = ms2$Intensity, 
                   Outpath = file.path(outPath, "SummedMS2.mgf"),
                   Title = paste("Scan numbers", paste(ScanNumbers, collapse = ", ")),
                   Pepmass = PrecursorsMzs,
                   Charge = PrecursorsCharges,
                   RT = "")
  
  result = isoforma_pipeline(Sequences = MultipleMods,
            SummedSpectra = ms2,
            MaxFragmentCharge = MaxFragmentCharge,
            ActivationMethod = "ETD",
            IonGroup = "c",
            CorrelationScore = CorrelationScore,
            PPMThreshold = PPMThreshold,
            IsotopeAlgorithm = "isopat")
  
  # Save results: -------
  
  if (!is.null(result[[1]]))
  {
    htmlwidgets::saveWidget(
      widget = annotated_spectrum_ptms_plot(SummedSpectra = ms2, IsoformaFragments = result[[1]]), #the plotly object
      file = file.path(outPath, "SpectrumAnnotated.html"), #the path & file name
      selfcontained = TRUE #creates a single html file
    )
    
    # Write annotated peaks:
    dat = NULL
    x = result[[1]]
    modNames = names(x)
    for(i in 1:length(x))
    {
      df = x[[i]]
      df$Modification = modNames[i]
      dat = rbind(dat, df)
    }
    write.csv(dat, file = file.path(outPath, "FragmentsPerModification.csv"), row.names = F, quote = F)
  }

  if (!is.null(result[[3]]))
    write.csv(result[[3]], file = file.path(outPath, "AbundanceMatrix.csv"), row.names = F, quote = F)
    
  if (!is.null(result[[4]]))
    write.csv(result[[4]], file = file.path(outPath, "Proportions.csv"), row.names = F, quote = F)
  
  if (!is.null(result[[5]]))
    ggplot2::ggsave(file.path(outPath, "Proportions.png"), result[[5]])
  
  
}


# ------------------------------
# Iterate rows in Targets.xslx file

cluster <- makeCluster(nCores)
registerDoParallel(cluster)

#for(k in 1:nrow(targets))
foreach(k = 1:nrow(targets),.packages=c("pspecterlib", "isoforma"),  .combine='c', .options.snow=list(preschedule=TRUE)) %dopar%
{
  mzml = targets$File[k]
  print(mzml)
  try(
  IsoFormaMain(mzMLFullPath = paste(mzMLFolder, mzml, sep = "/"),
               outPath = paste("Output/", sub(".mzml", "", mzml, ignore.case = TRUE), "_", targets$ProteoformGroup[k], sep = ""),
               UnmodSequence = targets$UnmodifiedSequence[k],
               Modification = targets$Modifications[k],
               RTStart = targets$RTStart[k],
               RTEnd = targets$RTEnd[k],
               parameters)
  )
}
stopCluster(cluster)
