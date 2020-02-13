# description: "Worm specific transformation file (part of src/transform.R)"
# author: "Manuel Belmadani"
# date: "29/04/2018"
# updated: "29/04/2018"

worm.tf.normalizeByWT  <- 
  # Descriptions
  function(X, variable, WT_VAR = "cePTEN(rf) + WT PTEN") {
    ## X=worm.df; variable="Worm.activity"
    scaleFunc <- function(Xi) {return( Xi / (right - left ) )}
    boundFunc <- function(Xi) {return (1 - Xi) }
    
    rownames(X) <- X$PTENVariant
    
    #left <- mean( X[ WT_CLASS , variable] )
    left = 0
    right <- mean( X[ WT_VAR, variable] )
    
    X[, variable] <- sapply(X[, variable], FUN = scaleFunc)
    #X[, variable] <- sapply(X[, variable], FUN = boundFunc)
    
    # # Sanity check:
    # stopifnot((((TRUE == cor(X$scaled, X$log2ratio)) == scaleFunc(right)) == (1 - scaleFunc(left))))
    
    return( X[, variable] )
  }


worm.tf.normalizeByCtrl <- 
  # Descriptions
  function(X, variable, WT_VAR = "Wildtype overexpression", EV_VAR = "cePTEN(rf)") {
    stop("UNTESTED")
    ## X=worm.df; variable="Worm.activity"
    scaleFunc <- function(Xi) {return( Xi / (right - left ) )}
    ### Omg this compiles boundFunc <- function(Xi) {return (Xi) - (right -1) }
    boundFunc <- function(Xi) {return (1 - Xi) }
    
    #left <- mean( X[ WT_CLASS , variable] )
    left = 0
    right <- mean( X[ EV_VAR, variable] )
    
    X[, variable] <- sapply(X[, variable], FUN = scaleFunc)
    X[, variable] <- sapply(X[, variable], FUN = boundFunc)
    
    # # Sanity check:
    # stopifnot((((TRUE == cor(X$scaled, X$log2ratio)) == scaleFunc(right)) == (1 - scaleFunc(left))))
    
    return( X[, variable] )
  }
