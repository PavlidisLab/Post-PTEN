# description: "Axon specific transformation file (part of src/transform.R)"
# author: "Manuel Belmadani"
# date: "11/05/2018"

axon.tf.normalizeByWT  <- 
  # Descriptions
  function(X, variable, WT_VAR = "WT", EV_VAR = "GFP") {
    ## X=worm.df; variable="Worm.activity"
    scaleFunc <- function(Xi) {return( Xi / (right - left ) )}
    boundFunc <- function(Xi) {return ( (right - left) - Xi ) }
    
    rownames(X) <- X$PTENVariant
    
    left <- 0
    right <- mean( X[ WT_VAR, variable] )
    
    X[, variable] <- sapply(X[, variable], FUN = scaleFunc)
    # X[, variable] <- sapply(X[, variable], FUN = boundFunc)
    
    return( X[, variable] )
  }
