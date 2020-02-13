## Instructions

Run `main.R` to process and the data and analysis.

## Content
_Top-level_  
  * README.md: You are here.
  * utils.R: General utilities (i/o, shorthands, contants)

_User directories_  
  * RawData: Data that's been provided as-is and is not ready for analysis in R (e.g. XLSX or DOCX files.)
  * Data: Processed data that has been preprocessed and ready for analysis.
  * Visualization: Utilities to render plots and other visualizations.
  
_Machine directories_  
  * res: Data that have been produced by scripts
  * src: Function definitions.
  * preprocessing: Scripts for loading files and applying correction (if any are required), and then saves or updates results to `res/`.
