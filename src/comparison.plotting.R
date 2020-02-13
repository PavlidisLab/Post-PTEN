from("src/plotting.R")
# description: "Plotting comparisons between multiple assays"
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "08/06/2018"

flyYeastScatterPlot <- 
  # Scatterplot of fly/yeast activity.
  # Example:
  # D=fly.normalized.df; Y=yeast.normalized.df; 
  function(Y,
           D,
           Y.metric="Yeast.activity",
           D.metric="Fly.activity",
           DIR="",
           showLabels=F){
    
    Combined <- merge(Y,
                      D[,colnames(D) != "Class"],
                      by="PTENVariant")
    
    Combined$Class <- as.factor(Combined$Class)
    
    
    colScale <- makeColScale(Combined, Variant="PTENVariant", type = "fill")  
    
    SORTED = order(rank(as.numeric(Combined[,D.metric]), ties.method = "first"))
    Combined = Combined[SORTED,]
    row.names(Combined) <- NULL
    
    PEARSON = signif(cor(Combined$Yeast.activity, Combined$Fly.activity, method = "pearson"),3)
    SPEARMAN = signif(cor(Combined$Yeast.activity, Combined$Fly.activity, method = "spearman"), 3)
    
    p <- ggplot(Combined, 
                aes(
                  x=Fly.activity,
                  y=Yeast.activity,
                  fill=factor(Combined[,"Class"])
                ))  +
      geom_point( size = DEFAULT_PT_SIZE*2, stat = "identity", shape=21, 
                  colour='black',
                  stroke=2) +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) +
      colScale + 
      TITLE_SIZES +
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      geom_abline(intercept = 0, slope = 1) +
      labs(title=paste0("Estimated variant activity (WT = 1)"), 
           subtitle = paste0("Pearson/Spearman correlation: ", PEARSON, " / ", SPEARMAN ),
           x= "Fly", 
           y = "Yeast" )
    
    
    
    if (showLabels){
      labelColours = classColors[  match(Combined$Class, names(classColors))  ]
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
                                  mapping=aes(label = as.character(Combined$PTENVariant)))
      
    }
    p
    # gggsave(paste0("-",metric), DIR=DIR)
    
    return(p)
  }
# flyYeastScatterPlot(D=fly.normalized.df, Y=yeast.normalized.df, showLabels = F)
# flyYeastScatterPlot(D=fly.normalized.df, Y=yeast.normalized.df, showLabels = T)

genericScatterPlot <- 
  # Scatterplot of fly/yeast activity.
  # Example:
  # G = CommonVariants; E=CommonVariants.err; G.metric1 = "mighell_score" ; G.metric2 = "Yeast.activity"; showLabels=F;
  function(G,
           G.metric1="mighell_score",
           G.metric2="Yeast.activity",
           E=NULL,
           DIR="",
           showLabels=F,
           showLines=c(0,1,0,1),
           PT_SIZE = DEFAULT_PT_SIZE,
           TITLE = "Estimated variant activity (WT = 1)",
           TITLE_SIZES.arg = TITLE_SIZES,
           SUBTITLE_PREFIX = "Pearson/Spearman correlation: ",
           ... ){
    
    if (!is.null(E)){
      
      for ( METRIC in c(G.metric1, G.metric2) ){
        if (METRIC %!in% colnames(E)){
          E <- cbind(E, PLACEHOLDER=as.numeric(G[,METRIC]))
          E$PLACEHOLDER[T] <- 0
          names(E)[ncol(E)] <- METRIC
          
        }
      }
      
      E <- E[!(is.na(G[,G.metric1]) | is.na(G[,G.metric2])),]
      E$Class <- factor(E$Class, levels = c(levels(E$Class), "Unclassified"))
      E$Class[is.na(E$Class)] <- "Unclassified"
      E$Class <- as.factor(E$Class)
      
    }
    
    G <- G[!(is.na(G[,G.metric1]) | is.na(G[,G.metric2])),]
    G$Class <- factor(G$Class, levels = unique(c(levels(G$Class), "Unclassified")) )
    G$Class[is.na(G$Class)] <- "Unclassified"
    G$Class <- as.factor(G$Class)
    colScale <- makeColScale(G, Variant="PTENVariant", type = "fill")
    colScaleColour <- makeColScale(G, Variant="PTENVariant", type = "colour")  
    
    
    # Set error bars
    XMAX = G[,G.metric1]
    XMIN = G[,G.metric1]
    YMAX = G[,G.metric2]
    YMIN = G[,G.metric2]
    if (!is.null(E)){
      XMAX[!is.na(E[,G.metric1])] = G[!is.na(E[,G.metric1]),G.metric1] + E[!is.na(E[,G.metric1]),G.metric1]
      XMIN[!is.na(E[,G.metric1])] = G[!is.na(E[,G.metric1]),G.metric1] - E[!is.na(E[,G.metric1]),G.metric1]
      
      YMAX[!is.na(E[,G.metric2])] = G[!is.na(E[,G.metric2]),G.metric2] + E[!is.na(E[,G.metric2]),G.metric2]
      YMIN[!is.na(E[,G.metric2])] = G[!is.na(E[,G.metric2]),G.metric2] - E[!is.na(E[,G.metric2]),G.metric2]
    }
    
    PEARSON = signif(cor(G[,G.metric1], G[,G.metric2], method = "pearson"),3)
    SPEARMAN = signif(cor(G[,G.metric1], G[,G.metric2], method = "spearman"), 3)
    
    p <- ggplot(data=G,
                mapping=aes(
                  x=G[,G.metric1],
                  y=G[,G.metric2],
                  
                  xmin=XMIN,
                  xmax=XMAX, 
                  
                  ymin=YMIN,
                  ymax=YMAX,
                  fill=factor(G[,"Class"]),
                  colour=factor(G[,"Class"])))  +
      geom_errorbar(width=0, size=1) +
      geom_errorbarh(size=1)  + 
      geom_point(
        size = PT_SIZE * 2,
        stat = "identity",
        shape=21,
        ...) + 
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) +
      colScale + colScaleColour + 
      TITLE_SIZES.arg + # geom_abline(intercept = 0, slope = 1) +
      labs(title=TITLE, 
           subtitle = paste0(SUBTITLE_PREFIX, PEARSON, " / ", SPEARMAN ),
           x = as.label(G.metric1), 
           y = as.label(G.metric2) ) 
    
    if (!is.null(E)){
      E[is.na(E[,G.metric1]), G.metric1] <- 2.0
      E[is.na(E[,G.metric2]), G.metric2] <- 2.0
    }
    
    if ( !isTRUE(showLines) ){
      p =  p + 
        geom_hline(yintercept = showLines[3], linetype="dashed") +
        geom_hline(yintercept = showLines[4], linetype="dashed") +
        geom_vline(xintercept = showLines[1], linetype="dashed") +
        geom_vline(xintercept = showLines[2], linetype="dashed")
    }
    
    
    if ( isTRUE(showLabels) ){
      labelColours = classColors[  match(G$Class, names(classColors))  ]
      textColours <- labelColours
      textColours[T] <- "white"
      textColours[labelColours %in% HIGH_CONTRAST] <- "black"
      
      p <- p + 
        ggrepel::geom_label_repel(data=G,
                                  mapping=aes(
                                    x=G[,G.metric1],
                                    y=G[,G.metric2],
                                    fill=factor(G[,"Class"]),
                                    label = as.character(G$PTENVariant)
                                  ),
                                  colour = textColours,
                                  segment.colour = "grey50",
                                  fontface="bold",
                                  box.padding = 0.35, 
                                  point.padding = 0.5,
                                  fill=labelColours)
      
    }
    # gggsave(paste0("-",metric), DIR=DIR)
    
    return(p)
  }

genericSinaPlot <- 
  # DESCRIPTION
  # Examples: D=yeast.normalized.df; 
  function(D, D.metric, groupBy="Class"){
    
    # Force to factor if not already.
    D[,groupBy] <- as.factor(D[,groupBy])
    
    valid_classes = table(D[,groupBy]) > 1
    valid_classes.names = names(valid_classes[valid_classes == T])
    
    D = D[D$Class %in% valid_classes.names,]
    n_groups <- length(valid_classes.names)    
    
    Y.TICKS=c(-0.5, 0.0, 0.5, 1.0, 1.5)
    par(mar = c(10,3,5,1) + 0.1)
    sinaplot(as.formula( paste(D.metric, "~", groupBy) ), 
             data = D,
             method = "count",
             pch = 20,
             xaxt = "n",
             yaxt = "n",
             col = as.vector(classColors[valid_classes.names]),
             labels=F,
             ann = FALSE, 
             bty = "n",
             cex=DEFAULT_PT_SIZE/2,
             ylim = c(-0.5, 1.7))
    axis(1, at = 1:n_groups, labels = FALSE) # X axis
    axis(2, at = Y.TICKS, labels=Y.TICKS, cex.axis=1.2) # Y Axis
    text(x = 1:n_groups,
         y = par()$usr[3] - 0.1 * (par()$usr[4] - par()$usr[3]),
         labels = valid_classes.names, srt = 25, xpd = TRUE, adj = 1,
         cex = 1.2)
  }

genericForestPlot <-
  # Comment
    # Example: M=combined ; E=combined.ci ;
    function(M,
             E,
             variant = "C124S",
             metrics = c("Axon.activity", "Fly.activity", "PSD.95.Density" , "Total.Dendrite.Length", "Worm.activity",  "Yeast.activity"),
             showSummary=F,
             sortByMean=F,
             Group = "PTENVariant",
             Colour = NULL,
             ...) {
      
      # Add Summary at the bottom.
      shapesVector = rep( 16, length(metrics))
      if (showSummary) {
        metrics <- c(metrics, "Summary")
        shapesVector = rep(16, length(metrics))
        shapesVector[length(shapesVector)] = 18
      }
      
      if (is.null(Colour)){
        variant.Class = ClassInformation[ClassInformation$PTENVariant == variant,]$Class
        Colour = classColors[as.character(variant.Class)]
        if (is.null(Colour) | is.na(Colour)){
            print("Defaulting to black...")
            Color = 'black'
          }
      }
      
      M.Genotypes = M[,Group]
      E.Genotypes = E[,Group]
      mean.effect  <- as.numeric(M[M.Genotypes == variant, metrics])
      err.lower <- as.numeric(mean.effect - E[E.Genotypes == variant, metrics])
      err.upper <- as.numeric(mean.effect + E[E.Genotypes == variant, metrics])

      
      df <- data.frame(metrics, mean.effect, err.lower, err.upper, shapesVector)
      
      if (sortByMean) {
        df$metrics <- factor(df$metrics, levels=rev(as.character(df$metrics[order(df$mean.effect)]) ))
        df$metrics <- lastlevel(df$metrics, "Summary") # Summary still gets last spot.
      } else { # Sort by class
        df$metrics <- factor(df$metrics, levels=rev(df$metrics) )
      } 
      
      fp <- ggplot(data=df, 
                   aes(x=as.label(df$metrics), y=mean.effect, ymin=err.lower, ymax=err.upper)
                   ) +
        geom_pointrange( size = 1, color=Colour, pch=shapesVector ) + 
        geom_hline(yintercept=0.5, lty=2) +  # add a dotted line at x=yintercept after flip
        coord_flip() +  # flip coordinates (puts labels on y axis)
        xlab("Assay") + ylab("Mean (95% CI)") + ggtitle(variant) +
        theme_bw()  # use a white background
      return(fp)  
    }
