if ( !exists("RESET_UTILS") ){
  RESET_UTILS = TRUE
}  

if ( RESET_UTILS == TRUE ) {
  
  # Options
  options(contrasts = rep ("contr.treatment", 2))
  options(error = utils::recover) # Proper error recovery/handling.
  # options(error = NULL) # Turn off debug
  
  # Utilities
  source("src/system.R")
  source("src/math.R")
  source("src/libraries.R")
  source("src/constants.R")
  source("src/io.R")
  source("src/sugar.R")
  source("src/plotting.R")
  source("src/transform.R")
}
RESET_UTILS = FALSE
