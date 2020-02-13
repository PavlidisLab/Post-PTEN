from("src/plotting.R")

PREFERED_ORDER.interactive = c("PTENVariant",
                               rev(c(
                                 "VAC7","VPS38","VPS30","VAM3","YPT7","VAM7","FIG4","VAC14",
                                 "Mighell2018",
                                 "Worm",
                                 "Rat.PSD95",
                                 "Rat.TDL",
                                 "Rat.Gephyrin",
                                 "Rat.Soma",
                                 "WB.Quant", 
                                 "Fowler.HTS", 
                                 "Fowler.IS"))
)



plotMOFA <- function(M, DIR="", PREFFERED_ORDER = NULL){
  
  #PREFERED_ORDER = colnames(M)[-1] #c( "Yeast", "Worm", "Rat")
  # PREFERED_ORDER = c("PTENVariant",
  #                    
  #                    #"△Yeast",
  #                    "Mighell2018",
  #                    "VAC14","VPS38","VPS30","VAM3","YPT7","VAM7","FIG4","VAC7",
  #                    "Worm",
  #                    "Rat.PSD95",
  #                    "Rat.TDL",
  #                    "Rat.Gephyrin",
  #                    "Rat.Soma",
  #                    #"WB.StdDev", 
  #                    "WB.Quant", 
  #                    "Fowler.HTS", 
  #                    "Fowler.IS")
  #               
  if (is.null(PREFFERED_ORDER)){
    PREFERED_ORDER = c("PTENVariant",rev(PREFERED_ORDER.interactive[-1]))
  } else if (is.logical(PREFFERED_ORDER) & PREFFERED_ORDER == F)  {
    PREFERED_ORDER = colnames(M)
  }
  
  M <- M[,PREFERED_ORDER]                     
  MEANORG <- rowMeans((M[ , PREFERED_ORDER[-1]]), na.rm = T)
  MEANORG.order <- order(rowMeans((M[, PREFERED_ORDER[-1]]), na.rm = T ))
  
  crosscompare.df.wMEANORG <- t(cbind(M, MOFA=MEANORG))
  colnames(crosscompare.df.wMEANORG) <- crosscompare.df.wMEANORG[1,]
  crosscompare.df.wMEANORG <- crosscompare.df.wMEANORG[-1,]
  crosscompare.df.wMEANORG <- apply(crosscompare.df.wMEANORG, MARGIN = 2, FUN = as.numeric)
  row.names(crosscompare.df.wMEANORG) <- c(PREFERED_ORDER[-1], "Mean")
  crosscompare.df.wMEANORG <- crosscompare.df.wMEANORG[rownames(crosscompare.df.wMEANORG) != "Mean",]

  printHeatmap <- function(){
    MATRIX <- crosscompare.df.wMEANORG
    # MATRIX["Rat.PSD95",] = MATRIX["Rat.PSD95",]
    # MATRIX["Rat.Gephyrin",] = MATRIX["Rat.Gephyrin",]
    MATRIX <- t(apply(MATRIX, MARGIN = 1, FUN = function(x){
      y <- x - min(x, na.rm = T)
      z <- y / max(y, na.rm = T)
    }))
    MATRIX = apply(MATRIX, MARGIN = 2, FUN = function(COL){
        COL - MATRIX[,"WT"] + 1.0
      })
            
            
    #ORDER = order(as.numeric(MATRIX["Yeast",]))
    cor(t(MATRIX), method = "pearson")
    
    heatmap.2(MATRIX[,] , 
              main = "Cross-assay comparison", # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins = c(8,12),     # widens margins around plot
              col=brewer.pal(7, "PRGn"),       # use on color palette defined earlier
              #col=my_palette,       # use on color palette defined earlier
              #breaks=col_breaks,    # enable color transition at specified limits
              # dendrogram="row",     # only draw a row dendrogram
              dendrogram='none',
              Colv = F, #ORDER,
              Rowv=F, 
              cexCol=2,
              cexRow = 2,
              na.color = 'darkgrey', 
              labCol = colnames(MATRIX),
              keysize=1
              )        
  }
  if (ENABLE_GGGSAVE){
    # Save to file.
    CairoSVG(file = paste0(DIR,"mofa-heatmap.svg"), width = 11, height = 7)
    printHeatmap()
    dev.off()
    
    png(file = paste0(DIR,"mofa-heatmap.png"),  width = 3, height=1)
    printHeatmap()
    dev.off()
  }
  else {
    # Just render it.
    printHeatmap()  
  }
} 
# plotMOFA(combined.numeric, PREFFERED_ORDER = F)
# plotMOFA(common.df)

yeast.pl.scatterplot <- 
  # Scatter plot to compare
  function(comparison.df, SENTINEL = "VAC14", showLabels = F, show.error=F){
    comparison.df$yeast = comparison.df[,p(SENTINEL,".Mean.normalized")]
    comparison.df$sd.yeast = comparison.df[,p(SENTINEL,".SD.normalized")]
    
    labelColours = classColors[  match(comparison.df$class, names(classColors))  ]
    p <- ggplot(data = comparison.df,
           mapping = aes(x = yeast,
                         y = mighell,
                         colour = class)) +
      colScale +
      geom_hline(yintercept = 0.5, linetype="dashed") +
      geom_vline(xintercept = 0.5, linetype="dashed")
    if (show.error){
      ERROR_SIZE = 10 # Reduce the factor by this magnitude.
      p <- p +
        geom_errorbar(mapping=aes(ymin=comparison.df$mighell - comparison.df$sem.mighell/ERROR_SIZE,
                                  ymax=comparison.df$mighell + comparison.df$sem.mighell/ERROR_SIZE), width=0, size=1, color='black') +
        
        geom_errorbarh(mapping=aes(xmin=comparison.df$yeast - comparison.df$sd.yeast/ERROR_SIZE,
                                   xmax=comparison.df$yeast + comparison.df$sd.yeast/ERROR_SIZE), size=1, color=labelColours )
    }

    p <- p + geom_point(mapping=aes(shape=comparison.df$isImputed), size=6)
    
    
    if (showLabels){
      labelColours = classColors[  match(comparison.df$class, names(classColors))  ]
      textColours <- labelColours
      textColours[T] <- "white"
      textColours[labelColours %in% HIGH_CONTRAST] <- "black"
      p <- p + 
        ggrepel::geom_label_repel(colour = textColours,
                                  segment.colour = "grey50",
                                  fontface="bold",
                                  box.padding = 0.35, 
                                  point.padding = 0.5,
                                  fill=labelColours,
                                  mapping=aes(label = as.character(comparison.df$variant)))
      
    } else {
      p <- p + theme_bw()  
    }
    THEME = theme(axis.title=element_text(size=20))
    return(p + THEME)
  }

########################################

yeastPointPlot  <- 
  # Example:
  # Y=yeast.normalized.df; Y.err=yeast.err.normalized.df; metric = "Yeast.activity";
  function(Y , Y.err = NULL, metric="Yeast.activity", DIR="", CI=F){
    
    ERR.factor = 1.0
    if (CI){
      ERR.factor = CI_FACTOR
    }
    
    dummy.gfp = c(0.0)
    dummy.err = c(0.0)
    Y <- rbind(Y, dummy.gfp)
    Y.err <- rbind(Y.err, dummy.gfp)
    
    # Relelvel variant names
    levels(Y$PTENVariant) <- c(levels((Y$PTENVariant)),  "Empty vector", "WT")
    levels(Y.err$PTENVariant) <- c(levels((Y.err$PTENVariant)),  "Empty vector", "WT")
    # Relevel classes
    levels(Y$Class) <- c( levels((Y$Class)),  "Control", "Wildtype overexpression")
    levels(Y.err$Class) <- c(levels((Y.err$Class)),  "Control", "Wildtype overexpression")
    
    Y[nrow(Y), c("PTENVariant", "Class")] <- c("Empty vector", "Control")

    Y.err$Class <- as.factor(Y.err$Class)
    Y$Class <- as.factor(Y$Class)
    
    colScale <- makeColScale(Y, Variant="PTENVariant")  
    if (is.numeric(metric)){
      metric <- colnames(Y)[metric]  
    }
    
    # Set minimum to 0.0
    SORTED = order(rank(as.numeric(Y[,metric]), ties.method = "first"))
    Y = Y[SORTED,]
    row.names(Y) <- NULL
    
    p <- ggplot(Y, 
                aes(
                  x=as.factor(1:nrow(Y)),
                  y=as.numeric(Y[,metric]),
                  colour=factor(Y[,"Class"])
                ))  +
      geom_point( size = DEFAULT_PT_SIZE*2, stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) +
      colScale + 
      TITLE_SIZES +
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_discrete(labels=as.character(Y$PTENVariant)) +
      labs(title=paste0("Yeast estimated variant activity (",metric,")"), 
           x= "Variants", 
           y = "Estimated activity (WT = 1)" )
    
    if (!is.null(Y.err)){
      # Set minimum to 0.0
      Y.err = Y.err[SORTED,]
      row.names(Y) <- NULL
      
      pd <- position_dodge(0.1) # move them .05 to the left and right
      
      p <-  p + geom_errorbar(mapping=aes(x=factor(rank( Y[,metric], ties.method = "first" )),
                                          ymin=Y[,metric] - (Y.err[,metric]*ERR.factor), 
                                          ymax=Y[,metric] + (Y.err[,metric]*ERR.factor), 
                                          colour=factor(Class)), 
                              width=DEFAULT_LN_WIDTH, 
                              size=DEFAULT_LN_SIZE, 
                              position=pd) 
    }
    
    # gggsave(paste0("wormPointPlot-",metric), DIR=DIR)
    p
    return(p)
  }

