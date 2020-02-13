# description: "Project wide constants"
# more: "Use for plotting or logical constants; in any case different set of constants should be swappable. E.g. a constant file could be made to generate plots with finer type, etc."
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "03/06/2018"
# from("utils.R") # Cannot be sourced for obvious reasons.

from <- 
  # This function doesn't actually do  anything other than clarity of sourced order and lets the user click to the hyperlinked source file.
  # Todo: if possible print something like "CurrentScript is from X"
  function(fromSource) {
    currentSource = "Current script"
    message = NULL
    a = tryCatch({
      currentSource = rstudioapi::getSourceEditorContext()$path
    }, warning = function(w) {
      print("[WARN] Warning encountered while running from().")
    }, error = function(e) {
      # Most likely RStudio isn't being used here.
      currentSource = strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2]
    }, finally = {
      message = paste(currentSource, "is called from", fromSource)
    })
    
    print(message)
    return(invisible(message))
  }

DEBUG_ON <- function() options(error = utils::recover)
DEBUG_OFF <- function() options(error = NULL)
