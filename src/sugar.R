# description: "Syntax sugar and shortcuts"
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
from("utils.R")

p <- paste0
ps <- 
  # Paste for paths with / separator.
  function(...){
    paste0(collapse = "/", ...)  
  }

indexOf <- 
  # Get the numeric index in dataframe by matching header by regex or string if fixed=T.
  function(string, df, fixed=F){
    grep(string, names(df), fixed = fixed)
  }

'%!in%' <- function(x,y){!('%in%'(x,y)) } # dplyr 'not in' statement.


trim <- function (x) gsub("^\\s+|\\s+$", "", x) # returns string w/o leading or trailing whitespace
trim.leading <- function (x)  sub("^\\s+", "", x) # returns string w/o leading whitespace
trim.trailing <- function (x) sub("\\s+$", "", x)  # returns string w/o trailing whitespace

howmuch <- function(Logical) length(which(Logical)) # Counts a logical, e.g. howmuch( "Bob" == c("Alice", "Bob", "Clyde", "Bob") ) would return 2

typeof.col <- function(Dataframe) apply(Dataframe, MARGIN = 2, FUN=typeof)

aaOrder <- 
  # Return the AA postion of the variant as a numeric 
  function(X) {
    X[!grepl("^[A-Za-z]+", X) ] <- NA
    X[!grepl("[A-Za-z]+$", X) ] <- NA
    aaTrimLeft <- function(X) gsub("^[a-zA-Z]+", "", X)
    aaTrimRight <- function(X) gsub("[a-zA-Z]+$", "", X)
    return(as.numeric(aaTrimRight(aaTrimLeft(X))))
  }

missenseOnly <-
  # Return missense variants only by filtering out frameshift "fs", or early stops (X, *)
  function(df, Index="PTENVariant") {
    return( df[ !grepl(df[,Index], pattern = "fs", fixed = T, ignore.case = F) & 
                !grepl(df[,Index], pattern = "X", fixed = T, ignore.case = F) &
                !grepl(df[,Index], pattern = "*", fixed = T, ignore.case = F), ] )
      
  }

truncatingOnly <-
  # Return non-missense variants only by filtering out non-frameshift "fs", or early stops (X, *)
  function(df, Index="PTENVariant") {
    return( df[ grepl(df[,Index], pattern = "fs", fixed = T, ignore.case = F) | 
                  grepl(df[,Index], pattern = "X", fixed = T, ignore.case = F) |
                  grepl(df[,Index], pattern = "*", fixed = T, ignore.case = F), ] )
    
  }

lastlevel <-
  # Reassign level to be the last.
  function(f, last) {
    if (!is.factor(f)) stop("f must be a factor")
    orig_levels = levels(f)
    if (! last %in% orig_levels) stop("last must be a level of f")
    new_levels = c(setdiff(orig_levels, last), last)
    factor(f, levels = new_levels)
  }

lastlevel = 
  # Same as other lastlevel but this one uses more standard methods and has more error checking.
  function (f, last, ...) {
    if (!is.factor(f)) stop("f must be a factor")
    lev <- levels(f)
    if (length(last) != 1L) 
      stop("'last' must be of length one")
    if (is.character(last)) 
      last <- match(last, lev)
    if (is.na(last)) 
      stop("'last' must be an existing level")
    nlev <- length(lev)
    if (last < 1 || last > nlev) 
      stop(gettextf("last = %d must be in 1L:%d", last, nlev), 
           domain = NA)
    factor(f, levels = lev[c(last, seq_along(lev)[-last])])
  }

subsetForBatch <-
  # Select non Genotype factors
  # TODO: This should be renamed to antipattern or something.
  function(X, pattern="Genotype", ...) {
    if (is.null(ncol(X))) { # Vector
      filter.genotype = !grepl(x=names(X), pattern = pattern,  ...)
      return( X[filter.genotype] )  
      
    } else { # Matrix
      filter.genotype = !grepl(x=colnames(X), pattern = pattern, ...)
      return( X[,filter.genotype] )  
      
    }
  }


source("src/sugar.labelling.R")
