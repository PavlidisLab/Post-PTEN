# description: "Subset of I/O utilities specific to yeast analysis"
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
from("utils.R")

yeast.io.load_multiple_csv <- 
  ## Given a path to multiple CSV files of equivalent header, combine them into a single dataframe.
  function(path, extension="*.txt") { 
    files <- dir(path, pattern = extension, full.names = TRUE)
    
    readWithFilename <- function(X){
      dat <- read.csv(X, header=TRUE, strip.white=TRUE, sep="\t")
      dat$fileName <- tools::file_path_sans_ext(basename(X))
      dat
    }
    
    tables <- lapply(files, readWithFilename)
    df <- data.frame(do.call(rbind, tables))
    
    df <- yeast.io.cleanup(df)
    
    return(df)
  }


yeast.io.cleanup <- 
  ## Return a cleaned up version of duplicate strings.
  function(dataset){
    
    # Rename screen to remove redundant plate name.
    dataset$Screen <- gsub("PTENPlateA--", "", dataset$Screen)
    # Fetch variant name from fileName
    dataset$Variant <- gsub("_set-1", "", dataset$fileName)
    
    # Remove unused columns
    relevant_columns <- c(  "Screen", "Gene.Name", "Plate", "Row", "Column", "Ctrl.Mean", "Exp.Mean", "Size.Difference", "fileName", "Variant" )
    names(dataset)
    dataset <- dataset[ c(relevant_columns) ]
    
    # Remove rows with 0 values
    dataset <- dataset[ (dataset$Ctrl.Mean != 0 & dataset$Exp.Mean != 0),  ]
    
    # Split Exp and Ctrl rows into two.
    dataset <- yeast.io.add_logical_rows(dataset)
    
    return(dataset)
  }


yeast.io.add_logical_rows <- 
  # Add the conceptual logical row to the array (Excludes the control row)
  function(df){
    
    # Add ratio
    df$Ratio <- df$Exp.Mean / df$Ctrl.Mean
    
    # Add a Size and Type to the dataframe.
    df$Size <- -1
    df$Type <- "UNASSIGNED"
    
    dfCtrl <- df
    dfExp <- df
    
    dfCtrl$Size <- dfCtrl$Ctrl.Mean
    dfExp$Size <- dfCtrl$Exp.Mean
    
    dfCtrl$Type <- "Ctrl"
    dfExp$Type <- "Exp"
    
    dfCtrl$LogicalRow <- dfCtrl$Row * 2 - 1
    dfExp$LogicalRow <- dfExp$Row * 2
    
    return_df <- rbind(dfCtrl, dfExp)
    
    return( return_df )
  }


yeast.io.load_pten_metadata <- function(filename){
  # TODO: Make this a standard format
  # filename = "PTEN_variants.csv"
  pten_meta <- read.csv(filename, header = TRUE, sep = "\t")
  
  for (name in names(pten_meta)){
    pten_meta[name] <- lapply(pten_meta[name], trim)
  }
  pten_meta$SNAP2 <- as.numeric(pten_meta$SNAP2)
  pten_meta$CADD <- as.numeric(pten_meta$CADD)
  
  return(pten_meta)
}

yeast.io.match_pten_metadata <- 
  
  # TODO: Decouple data "fixing" steps from metadata annotation.
  # Match PTEN metadata from external source (yeast.io.load_pten_metadata())
  function(df, PTEN_METADATA){
  
  pten_meta <- yeast.io.load_pten_metadata(PTEN_METADATA)
  # head(df)
  # head(pten_meta)
  # length(unique(pten_meta$PTENVariant) )
  # length(unique(df$Variant)  )

  setdiff(unique(pten_meta$PTENVariant) ,unique(df$Variant))
  setdiff(unique(df$Variant),unique(pten_meta$PTENVariant))
  
  # Join metadata where mutation matches
  df <- merge(x = df, y = pten_meta, by.x="Variant", by.y = "PTENVariant", all.x=TRUE, sort = TRUE)
  
  
  VariantStatus <- vector(mode="character", length=nrow(df))
  df <- cbind(df, VariantStatus)
  df$VariantStatus <- df$Class
  
  
  df$VariantStatus[ df$Type == "Ctrl" ] <- "ControlSpot"
  df$VariantStatus[ df$Type != "Ctrl" ] <- df$Variant[ df$Type != "Ctrl" ]
  
  # add 'block' data
  df$blockRow = ceiling(df$LogicalRow/4)
  df$blockCol = ceiling(df$Column/4)
  df$block = factor(df$blockRow + 8*(df$blockCol - 1))
  df$edge = df$LogicalRow == 1 | df$Column == 1 | df$LogicalRow == 32 |  df$Column == 48
  
  
  # These we don't have metadata for.
  unique(df$Variant)[ unique(df$Variant) %!in% (pten_meta$Variant) ]
  
  # Fix classes
  df$Class[is.na(df$Class)] <- "Unclassified"
  df$Class[grep("^WT", df$Variant)] <- STR_WT
  df$Class[grep("pEGH", df$Variant)] <- STR_EV
  
  # Finally, fix types.
  df$SNAP2 <- as.numeric(df$SNAP2)
  df$CADD <- as.numeric(df$CADD)
  
  
  
  return(df)
  }

yeast.io.fixTypos <-
  # Fix known typos in the raw data
  function(x, batch=1){
    mistakes   <- c("A125P", "G285X", "G38E",  "N356Q", "Q346R", "R357S", "V343M", "C136fs",  "C136M", "T240X", "T68X")
    correction <- c("A126P", "E285X", "G36E",  "N356D", "Q396R", "P357S", "V343L", "I135fs", "I135fs", "Y240X", "Y68X")
    
    x[ !is.na(match(x, mistakes))] = correction[ match(x, mistakes)[ !is.na(match(x, mistakes))]  ]
    
    x[x %in% mistakes] <- correction[ mistakes %in% unique(x[x %in% mistakes]) ]
    
    return(x)
  }

yeast.io.cast <- 
  # Create an aggregated matrix of value "variable" using funcion "FUN".
  # Used to cast the data frame into mean rows by sentinel (Gene.Name)
  function(X, variable="log2ratio", FUN=mean) {
  
    b <- dcast(X, value.var=variable, fun.agg=FUN, Variant ~ Gene.Name) # Assuming this must be done on a Variant/Sentinel level
    bnames <- b[,1] # Names
    mtrx <- as.matrix(b[,-1])
    row.names(mtrx) <- bnames
    
    return(mtrx)
  }

yeast.io.compactAggregatedDF <- 
  # Squash the dataframe based on non-duplicate Sentinel, Variant, Mean.normalized
  function(X){
    return (X[!duplicated(X[,c("Gene.Name", "Variant", "Mean.log2ratio", "SD.log2ratio", "Mean.normalized", "SD.normalized", "Class")]), ]
            [,c("Gene.Name", "Variant", "Mean.log2ratio", "SD.log2ratio", "Mean.normalized", "SD.normalized", "Class", "SNAP2", "CADD", "Mean.zYeast",  "SD.zYeast")])
  }

