# description: "Transformation of data (scaling, normalizing, summary statistics.)"
# author: "Manuel Belmadani"
# date: "29/04/2018"
# updated: "29/04/2018"
from("utils.R")

###########################################################
## Core
###########################################################

normalize01 <- 
  # Normalize a vector to be strictly between 0 and 1.
  function(x){
    y <- x - min(x, na.rm = T)
    z <- y / max(y, na.rm = T)
    return(z)
  }

correctErrorBars <-
  # Apply correction on error bars based on an A -> B prior transformation;- 
  # rescales to the same magnitude
  # Examples:
  # E = fly.err.normalized.df$Fly.activity ; Before = fly.df$Fly.activity ;  After = fly.normalized.df$Fly.activity; Mask = colnames(Before)
  # E = worm.err.normalized.df; Before = worm.df; After = worm.normalized.df; Mask = "Worm.activity"
  function(E, Before, After, Mask=NULL) {
    E.correction <- E
    
    if (is.null(Mask)){
      # Use ALL columns
      Mask = colnames(Before) 
      
    } else if (is.vector(Mask)) {
      # Do nothing. Mask is already a vector
    }  else if (is.character(Mask)) {
      # Convert Mask to a vector
      Mask = c(Mask)
    } else {
      stop("Unknown Mask type in error correction transform!")  
    }
    
    sapply(Mask,
           function(x) {
             # Scale by the ratio of deltaB : deltaA
             scaling_factor = (max(After[,x])-min(After[,x])) / 
               (max(Before[,x])-min(Before[,x]))
             E.correction[,x] <<- E.correction[,x] * scaling_factor
           })
    return(E.correction)
  }

changeClass <-
  #
  function(DF, newclasses){
    DF <- merge(DF, newclasses, by="PTENVariant")
    DF$Class[!is.na(DF$Class.x)] = as.character(DF$Class.x[!is.na(DF$Class.x)])
    DF$Class[!is.na(DF$Class.y)] = as.character(DF$Class.y[!is.na(DF$Class.y)])
    DF$Class = as.factor(DF$Class)    
    DF$Class.x = NULL
    DF$Class.y = NULL
    
    return(DF)
  }

adjust.covariates <-
  function(x, model, pattern){
    
    coefs = coefficients(model)
    coefs[is.na(coefs)] = 0
    
    design = model.matrix(model)  
    
    return(
      t(x - subsetForBatch(coefs, pattern = pattern) %*% t(subsetForBatch(design, pattern = pattern)) )
    )
  }

adjust.data.for.batch <-
  # Batch correction using the model matrix; subtract out the variance attributatble to day/plate effects.
  # Defaults: variant.name="Genotype"; sd.name=NULL; adjs.name=NULL; adjs.sds.name=NULL
  function(x, model, value.name, variant.name="Genotype", pattern=NULL, sd.name=NULL, adjs.name=NULL, adjs.sds.name=NULL){
    
    if ( is.null(sd.name) ){
      sd.name = paste0("sd.", value.name)
    }
    sem.name = paste0("sem.", value.name)
    
    if ( is.null(adjs.name) ){ 
      adjs.name = paste0(value.name, ".adj")
    } else if ( !grepl(value.name, adjs.name) ) {
      stop("[ERROR] value.name must be asubstring of adjs.name.")
    }
    if ( is.null(adjs.sds.name) ){
      adjs.sds.name = paste0("sd.", adjs.name)
    }
    adjs.sem.name = paste0("sem.", adjs.name)
    
    if (is.null(pattern)){
      pattern = variant.name  
    }
    
    # Adjust for covariates
    x[, adjs.name] <- as.numeric(adjust.covariates(x[,value.name], model, pattern))
    
    filter.agg = grepl(x = names(x), pattern= value.name ) 
    
    computeCI <- function(s,n) qnorm(0.975)*s / sqrt(n)
    computeSEM <- function(s,n) s / sqrt(n)
    
    # aggregate(as.numeric(as.matrix(x[,"Measure"])),
    #           by = list(Condition = as.factor(x[,variant.name]) ),
    #           FUN=mean)
    
    mean.data = Filter(function(x) !all(is.na(x)),
                       aggregate(data.frame(x[,filter.agg]),
                          by = list(variant.name = x[,variant.name] ),
                          FUN=mean))
    names(mean.data)[1] = variant.name
    
    sd.data = Filter(function(x) !all(is.na(x)),
                     aggregate(x[,filter.agg], 
                        by = list(variant.name = x[,variant.name]),
                        FUN=sd))
    names(sd.data)[1] = variant.name
    names(sd.data) = paste0("sd.", names(sd.data))
    
    n.data = aggregate(x[,filter.agg], 
                       by = list(Genotype = x[,variant.name]),
                       FUN=length)
    names(n.data)[1] = variant.name
    names(n.data) = paste0("N.", names(n.data)) 
    
    mean.comparison.with.adjusted = merge(mean.data, 
                                          sd.data, 
                                          by.x = variant.name, 
                                          by.y = paste0("sd.", variant.name) )
    
    mean.comparison.with.adjusted = merge(mean.comparison.with.adjusted, 
                                          n.data, 
                                          by.x = variant.name, 
                                          by.y = paste0("N.", variant.name) )
    
    
    mean.comparison.with.adjusted[,sem.name] = computeSEM(mean.comparison.with.adjusted[,sd.name], mean.comparison.with.adjusted[,paste0("N.", value.name)])
    mean.comparison.with.adjusted[,adjs.sem.name] = computeSEM(mean.comparison.with.adjusted[,adjs.sds.name], mean.comparison.with.adjusted[,paste0("N.", value.name)])
    
    return(mean.comparison.with.adjusted)
  }

adjust.data.for.batch.lme4 <-
  # Batch correction using the model matrix; subtract out the variance attributatble to day/plate effects.
  # Defaults: variant.name="Genotype"; sd.name=NULL; adjs.name=NULL; adjs.sds.name=NULL
  function(x, model, model.rnef, value.name, rnef.names, variant.name="Genotype", filter.agg = NULL, pattern=NULL, sd.name=NULL, adjs.name=NULL, adjs.sds.name=NULL, agg.by.list=F){
    
    if ( is.null(sd.name) ){
      sd.name = paste0("sd.", value.name)
    }
    sem.name = paste0("sem.", value.name)
    
    if ( is.null(adjs.name) ){ 
      adjs.name = paste0(value.name, ".adj")
    } else if ( !grepl(value.name, adjs.name) ) {
      stop("[ERROR] value.name must be asubstring of adjs.name.")
    }
    if ( is.null(adjs.sds.name) ){
      adjs.sds.name = paste0("sd.", adjs.name)
    }
    adjs.sem.name = paste0("sem.", adjs.name)
    
    if (is.null(pattern)){
      pattern = variant.name  
    }
    
    # Adjust for covariates
    # model.rnef.df = data.frame(model.rnef)
    # x[, adjs.name] <- x[,value.name]
    # for ( effect in rnef.names ) {
    #   x[, adjs.name] <- x[,value.name] - model.rnef.df[model.rnef.df$grpvar == effect, "condval"]
    # }
    x[, adjs.name] <- getME(model, "X") %*% getME(model, "beta") + residuals(model)
    
    if (is.null(filter.agg)){
      filter.agg = grepl(x = names(x), pattern= value.name )  
    } else {
      filter.agg = grepl(x = names(x), pattern= filter.agg )  
    }
    
    computeCI <- function(s,n) qnorm(0.975)*s / sqrt(n)
    computeSEM <- function(s,n) s / sqrt(n)
    
    if (agg.by.list){
      agg.by = list(Genotype = x[,variant.name])
    } else {
      agg.by = x[,variant.name]
    }
    
    mean.data = aggregate(data.frame(x[,filter.agg]),
                          by = agg.by,
                          FUN=mean)
    
    
    names(mean.data)[1] = variant.name
    
    sd.data = aggregate(data.frame(x[,filter.agg]),
                        by = agg.by,
                        FUN=sd)
    
    names(sd.data)[1] = variant.name
    names(sd.data) = paste0("sd.", names(sd.data))
    
    n.data = aggregate(data.frame(x[,filter.agg]),
                       by = agg.by,
                       FUN=length)
    
    names(n.data)[1] = variant.name
    names(n.data) = paste0("N.", names(n.data)) 
    
    mean.comparison.with.adjusted = merge(mean.data, 
                                          sd.data, 
                                          by.x = variant.name, 
                                          by.y = paste0("sd.", variant.name) )
    
    mean.comparison.with.adjusted = merge(mean.comparison.with.adjusted, 
                                          n.data, 
                                          by.x = variant.name, 
                                          by.y = paste0("N.", variant.name) )
    
    
    mean.comparison.with.adjusted[,sem.name] = computeSEM(mean.comparison.with.adjusted[,sd.name], mean.comparison.with.adjusted[,paste0("N.", value.name)])
    mean.comparison.with.adjusted[,adjs.sem.name] = computeSEM(mean.comparison.with.adjusted[,adjs.sds.name], mean.comparison.with.adjusted[,paste0("N.", value.name)])
    
    return(mean.comparison.with.adjusted)
  }


add.scaled.01.means <- 
  #  df = adjusted.means; measurement = "Chemotaxis_Index.adj"; ctrl.null = PTEN_NULL; ctrl.wt = PTEN_WT;
  function(df, measurement, ctrl.null, ctrl.wt, select.by="Genotype") {
    original.means = df[,measurement]
    scaled.means = original.means - df[df[,select.by] == ctrl.null, measurement]
    scaled.means = scaled.means / scaled.means[df[,select.by] == ctrl.wt] 
    
    scaling_factor = (max(scaled.means)-min(scaled.means)) / 
      (max(original.means)-min(original.means))
    
    sd.name = paste0("sd.", measurement)
    sem.name = paste0("sem.", measurement) 
    
    scaled.sd = df[,sd.name] * scaling_factor
    scaled.sem = df[,sem.name] * scaling_factor
    df.scaled = data.frame(scaled.means,
                           scaled.sd,
                           scaled.sem)
    names(df.scaled) = gsub(pattern = "scaled", replacement = paste0(measurement, ".scaled"), x = names(df.scaled))
    
    cbind(df,df.scaled)
  }

###########################################################
###########################################################
source("src/yeast.transform.R")
source("src/fly.transform.R")
source("src/worm.transform.R")
source("src/rat.transform.R")
source("src/axon.transform.R")
