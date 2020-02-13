# description: "I/O utilities"
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
from("utils.R")

###########################################################
## Core
###########################################################
writeObject <- 
  # For a given dataframe, create a file using its variable name and store it in DATA.DIR
  function(Dataframe, DIR=DATA.DIR, filename=NULL, row.names = T, ...){
    if (is.null(filename)){
      filename <-  deparse(substitute(Dataframe))
    }
    filepath <- paste0(DIR, filename)
    write.table(Dataframe, file=filepath, sep="\t", row.names = row.names, ...)
  }

readMatrix <- 
  # Doesn't assume the path starts with "data/"
  function(path){
    return(as.matrix(read.table(path, sep = "\t")))
  }

readObjectAsMatrix <- 
  # Uses common DATA.DIR
  function(path, ...){
    return(as.matrix(read.table(paste0(DATA.DIR, path), sep = "\t", ...)))
  }

readObjectAsTable <- 
  # Uses common DATA.DIR
  function(path, ...){
    return(read.table(paste0(DATA.DIR,path), sep = "\t", ...))
  }

############################################################
############################################################
source("src/yeast.io.R")
