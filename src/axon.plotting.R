from("src/plotting.R")

axonPointPlot  <- 
  # Example:
  # A=axon.normalized.df; A.err=axon.err.normalized.df; metric = "Axon.activity";
  function(A , A.err = NULL, metric="Axon.activity", DIR="", CI=F){
    
    ERR.factor = 1.0
    if (CI){
      ERR.factor = CI_FACTOR
    }
    
    A <- rbind(A, c(NA, 0, NA))
    levels(A$PTENVariant) <- unique(c(as.character(A$PTENVariant),  "GFP"))
    A[nrow(A), c(1,3)] <- c("GFP", "Control")
    
    A$Class <- factor(A$Class)
    colScale <- makeColScale(A, Variant="PTENVariant")  
    if (is.numeric(metric)){
      metric <- colnames(A)[metric]  
    }
    
    # Set minimum to 0.0
    SORTED = order(rank(as.numeric(A[,metric]), ties.method = "first"))
    A = A[SORTED,]
    row.names(A) <- NULL
    
    p <- ggplot(A, 
                aes(
                  x=as.factor(1:nrow(A)),
                  y=as.numeric(A[,metric]),
                  colour=factor(A[,"Class"])
                ))  +
      geom_point( size = DEFAULT_PT_SIZE*2, stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=14)) +
      colScale + 
      TITLE_SIZES +
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_discrete(labels=as.character(A$PTENVariant)) +
      labs(title=paste0("Axonal outgrowth estimated variant activity (",metric,")"), 
           x= "Variants", 
           y = "Estimated activity (WT = 1)" )
    
    if (!is.null(A.err)){
      A.err <- rbind(A.err, c(NA, 0, NA)) # TODO find a better may to pass sd/error
      A.err = A.err[SORTED,]
      A.err$Class <- as.factor(A$Class) # Transfer class
      levels(A.err$PTENVariant) <- A$PTENVariant
      
      pd <- position_dodge(0.1) # move them .05 to the left and right
      
      p <-  p + geom_errorbar(mapping=aes(x=factor(rank( A[,metric], ties.method = "first" )),
                                          ymin=A[,metric] - (A.err[,metric]*ERR.factor), 
                                          ymax=A[,metric] + (A.err[,metric]*ERR.factor), 
                                          colour=factor(Class)), 
                              width=DEFAULT_LN_WIDTH, 
                              size=DEFAULT_LN_SIZE, 
                              position=pd) 
    }
    
    # gggsave(paste0("wormPointPlot-",metric), DIR=DIR)
    p
    return(p)
  }
