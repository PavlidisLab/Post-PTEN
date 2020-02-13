
wormPointPlot <- 
  # Example:
  # W=worm.normalized.df; W.err=worm.err.normalized.df; metric = "Worm.activity"
  function(W , W.err = NULL, metric="Worm.activity", DIR="", CI=F){
    
    ERR.factor = 1.0
    if (CI){
      ERR.factor = CI_FACTOR
    }
    
    #W = cbind(PTENVariant = rownames(W), W)
    # W <- rbind(W, `cePTEN(rf)`=c(0,  0, 0, NA))
    W["cePTEN(rf)", "Class"] <- "Control"
    # W = cbind(PTENVariant=row.names(W),
    #           W)
    colScale <- makeColScale(W, Variant="PTENVariant")  
    if (is.numeric(metric)){
      metric <- colnames(W)[metric]  
    }
    
    # Set minimum to 0.0
    SORTED = order(rank(as.numeric(W[,metric]), ties.method = "first"))
    W = W[SORTED,]
    row.names(W) <- NULL
    
    p <- ggplot(W, 
                aes(
                  x=as.factor(1:nrow(W)),
                  y=as.numeric(W[,metric]),
                  colour=factor(W[,"Class"])
                ))  +
      geom_point( size = DEFAULT_PT_SIZE*2, stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) +
      colScale + 
      TITLE_SIZES +
      # geom_hline(yintercept = W[W$Class == STR_WT,metric] - W.err[W$Class == STR_WT,metric], linetype="dashed") +
      # geom_hline(yintercept = W[W$Class == STR_WT,metric] + W.err[W$Class == STR_WT,metric], linetype="dashed") +
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      
      scale_x_discrete(labels=as.character(W$PTENVariant)) +
      labs(title=paste0("Worm estimated variant activity (",metric,")"), 
           x= "Variants", 
           y = "Estimated activity (N2 WT = 1)" )
    
    if (!is.null(W.err)){
      W.err <- rbind(W.err, `cePTEN(rf)`=c(0,  0, 0, NA))
      W.err[nrow(W.err), 4] <- "Control"
      W.err = cbind(PTENVariant=row.names(W.err),
                W.err)
      W.err = W.err[SORTED,]
      
      pd <- position_dodge(0.1) # move them .05 to the left and right
      
      p <-  p + geom_errorbar(mapping=aes(x=factor(rank( W[,metric], ties.method = "first" )),
                                          ymin=W[,metric] - (W.err[,metric] * ERR.factor ), 
                                          ymax=W[,metric] + (W.err[,metric] * ERR.factor ), 
                                          colour=factor(Class)), 
                              width=DEFAULT_LN_WIDTH, 
                              size=DEFAULT_LN_SIZE, 
                              position=pd) 
    }
    
    # gggsave(paste0("wormPointPlot-",metric), DIR=DIR)
    p
    return(p)
}
