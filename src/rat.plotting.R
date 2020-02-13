from("utils.R")

ratPointPlot <- 
  # Example:
  # R=rat.normalized.df; R.err=rat.err.normalized.df; metric="PSD.95.Density";
  function(R , R.err = NULL, metric="PSD.95.Density", DIR="", CI=F){
    
    ERR.factor = 1.0
    if (CI){
      ERR.factor = CI_FACTOR
    }
    
    R <- rbind(R, `GFP Control`=c(NA, 0,  0, 0, 0, NA))
    R[nrow(R), c("PTENVariant", "Class")] <- c("GFP Control", "Control")
    
    colScale <- makeColScale(R, Variant="PTENVariant")  
    if (is.numeric(metric)){
      metric <- colnames(R)[metric]  
    }
    
    
    # Set minimum to 0.0
    SORTED = order(rank(as.numeric(R[,metric]), ties.method = "first"))
    R = R[SORTED,]
    row.names(R) <- NULL
    
      ############ STILL SCREWED UP
    p <- ggplot(R, 
                aes(
                  x=as.factor(1:nrow(R)),
                  y=as.numeric(R[,metric]),
                  colour=factor(R[,"Class"])
                    ))  +
      geom_point( size = DEFAULT_PT_SIZE*2, stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) +
      colScale + 
      TITLE_SIZES +
      # geom_hline(yintercept = R[R$Class == STR_WT,metric] - R.err[R$Class == STR_WT,metric], linetype="dashed") +
      # geom_hline(yintercept = R[R$Class == STR_WT,metric] + R.err[R$Class == STR_WT,metric], linetype="dashed") +
      
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_discrete(labels=as.character(R$PTENVariant)) +
      labs(title=paste0("Rat estimated variant activity (",metric,")"), 
           x= "Variants", 
           y = "Estimated activity" )
    
    if (!is.null(R.err)){
      pd <- position_dodge(0.1) # move them .05 to the left and right
      
      R.err <- rbind(R.err, `GFP Control`=c(NA, 0,  0, 0, 0, NA))
      R.err[nrow(R.err), c("PTENVariant", "Class")] <- c("GFP Control", "Control")
      R.err = R.err[SORTED,]
      
      p <-  p + geom_errorbar(mapping=aes(x=factor(rank( R[,metric], ties.method = "first" )),
                                          ymin=R[,metric] - (R.err[,metric] * ERR.factor), 
                                          ymax=R[,metric] + (R.err[,metric] * ERR.factor), 
                                          colour=factor(Class)), 
                              
                              width=DEFAULT_LN_WIDTH, 
                              size=DEFAULT_LN_SIZE, 
                              position=pd) 
    }
    
    #gggsave(paste0("ratPointPlot-",metric), DIR=DIR)
    p
    return(p)
  }
# ratPointPlot(rat.normalized.df)
