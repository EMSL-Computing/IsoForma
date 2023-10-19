
# IsoForma demo
An R script implemented as a demo for processing multiple LC-MS/MS files accompanying the manuscript: "IsoForma: An R Package for Quantifying and Visualizing Positional Isomers in Top-Down LC-MS/MS Data".
The source code of the IsoForma package and documentation in vignettes are available at https://github.com/EMSL-Computing/isoforma-lib.

### USAGE (Windows)

1. Download and install R: https://ftp.osuosl.org/pub/cran/
2. Download, install, and open RStudio (Free version): https://www.rstudio.com/products/rstudio/download/
3. If re-installing, remove any previous versions of pspecter
	remove.packages("pspecterlib")
4. Clone/download IsoForma-paper repository
5. Open IsoForma_demo.R in RStudio
6. Install packages executing lines 6-10
7. Update the path in line 19 according to the location of the data in your computer
8. Execute all code in the script. Press Ctrl+Alt to select all code and Ctrl+Enter to execute it.

### Example data

* Parameters.xlsx
* Targets.xlsx
* RawData folder contains 2 mzML files converted using MSConvert (v. 3.0.22139-aa8be89) using parameters: scan time 5400-7800 seconds, vendor peak picking and ms levels 1-2.

### Reference

If you use IsoForma or any portions of this code please cite: Degnan et al. "IsoForma: An R Package for Quantifying and Visualizing Positional Isomers in Top-Down LC-MS/MS Data". Submitted to Journal of Proteome Research.
