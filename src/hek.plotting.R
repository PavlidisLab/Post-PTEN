from("src/plotting.R")


hekPointPlot  <- 
  # Example:
  # H=protein.df; H.err=protein.err; metric = "Hek.stability";
  function(H , H.err = NULL, metric="Hek.stability", PT_SIZE=DEFAULT_PT_SIZE*2, DIR="", CI=F){
    
    ERR.factor = 1.0
    if (CI){
      ERR.factor = CI_FACTOR
    }
    
    dummy.gfp = c(0.0)
    dummy.wt = c(1.0)
    dummer.err = c(0.0)
    # length(dummy.wt) <- ncol(H)
    # length(dummy.gfp) <- ncol(H)
    
    H <- rbind(H, dummy.gfp)
    H <- rbind(H, dummy.wt)
    
    # Relelvel variant names
    levels(H$PTENVariant) <- c(levels((H$PTENVariant)),  "GFP", "WT")
    levels(H.err$PTENVariant) <- c(levels((H.err$PTENVariant)),  "GFP", "WT")
    # Relevel classes
    levels(H$Class) <- c( levels((H$Class)),  "Control", "Wildtype overexpression")
    levels(H.err$Class) <- c(levels((H.err$Class)),  "Control", "Wildtype overexpression")
    
    
    H[nrow(H)-1, c("PTENVariant", "Class")] <- c("GFP", "Control")
    H[nrow(H),   c("PTENVariant", "Class")] <- c("WT", "Wildtype overexpression")
    
    H.err$Class <- as.factor(H.err$Class)
    H$Class <- as.factor(H$Class)
    
    colScale <- makeColScale(H, Variant="PTENVariant")  
    if (is.numeric(metric)){
      metric <- colnames(H)[metric]  
    }
    
    # Set minimum to 0.0
    SORTED = order(rank(as.numeric(H[,metric]), ties.method = "first"))
    H = H[SORTED,]
    row.names(H) <- NULL
    
    p <- ggplot(H, 
                aes(
                  x=as.factor(1:nrow(H)),
                  y=as.numeric(H[,metric]),
                  colour=factor(H[,"Class"])
                ))  +
      geom_point( size = PT_SIZE, stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) +
      colScale + 
      TITLE_SIZES +
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_discrete(labels=as.character(H$PTENVariant)) +
      labs(title=paste0("PTEN Hek stability assay"), 
           x= "Variants", 
           y = "Estimated activity (WT = 1)" )
    
    if (!is.null(H.err)){
      H.err <- rbind(H.err, c(NA, 0, NA)) # TODO find a better may to pass sd/error
      H.err = H.err[SORTED,]
      
      pd <- position_dodge(0.1) # move them .05 to the left and right
      
      p <-  p + geom_errorbar(mapping=aes(x=factor(rank( H[,metric], ties.method = "first" )),
                                          ymin=H[,metric] - (H.err[,metric]*ERR.factor), 
                                          ymax=H[,metric] + (H.err[,metric]*ERR.factor), 
                                          colour=factor(Class)), 
                              width=DEFAULT_LN_WIDTH, 
                              size=DEFAULT_LN_SIZE, 
                              position=pd) 
    }
    
    # gggsave(paste0("wormPointPlot-",metric), DIR=DIR)
    p
    return(p)
  }
  
