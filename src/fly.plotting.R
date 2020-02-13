from("src/plotting.R")

flyPointPlot <- 
  # Example:
  # D=fly.normalized.df; D.err=fly.err.normalized.df; metric = "Fly.activity";
  function(D , D.err = NULL, metric="Fly.activity", DIR="", CI=F, LETHAL=c()){
    
    ERR.factor = 1.0
    if (CI){
      ERR.factor = CI_FACTOR
    }
    
    ## Manually plot variants in the LETHAL vector.
    D <- D[  D$PTENVariant %!in% LETHAL, ]
    D.err <- D.err[  D.err$PTENVariant %!in% LETHAL, ]
    LETHAL_PENALTY = max(D[,metric]) + 0.1 # Top of the plot
    for (lethal_variant in LETHAL){
      LETHAL_ROW = data.frame(PTENVariant = lethal_variant,
                               Fly.activity = LETHAL_PENALTY,
                               Class = ClassInformation$Class[ClassInformation$PTENVariant==lethal_variant])
      D <- rbind(D, LETHAL_ROW)
      
      # If the standard error 
      # Question: Should this just be forced? e.g. Lethal -> no std.err
      
      if (
        !is.null(D.err) & # If error vector provided
        length(intersect(c(lethal_variant),D.err$PTENVariant)) == 0 # and variant NOT in error vector
          ){
        LETHAL_ROW_ERROR = data.frame(
          PTENVariant=lethal_variant, 
          Fly.activity=0.0, 
          Class = ClassInformation$Class[ClassInformation$PTENVariant == lethal_variant]
        )
        D.err <- rbind(D.err, LETHAL_ROW_ERROR)
      }
      
    }
    
    if (F & "attp2" %!in% D$PTENVariant){
      # TODO: Dead code
      D <- rbind(D, c(NA, 0,  0, 0, NA))
      levels(D$PTENVariant) <- unique(c(as.character(D$PTENVariant),  "attp2"))
      levels(D.err$PTENVariant) <- unique(c(as.character(D.err$PTENVariant),  "attp2"))
      D[nrow(D), c(1,3)] <- c("attp2", "Control")
      D.err$Class <- as.factor(D.err$Class)
      
    }
    D$Class <- as.factor(D$Class)
    
    
    colScale <- makeColScale(D, Variant="PTENVariant")  
    if (is.numeric(metric)){
      metric <- colnames(D)[metric]  
    }
    
    # Set minimum to 0.0
    SORTED = order(rank(as.numeric(D[,metric]), ties.method = "first"))
    D = D[SORTED,]
    row.names(D) <- NULL
    
    p <- ggplot(D, 
                aes(
                  x=as.factor(1:nrow(D)),
                  y=as.numeric(D[,metric]),
                  colour=factor(D[,"Class"])
                ))  +
      geom_point( size = DEFAULT_PT_SIZE*2, stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=10)) +
      colScale + 
      TITLE_SIZES +
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      scale_x_discrete(labels=as.character(D$PTENVariant)) +
      labs(title=paste0("Fly estimated variant activity (",metric,")"), 
           x= "Variants", 
           y = "Estimated activity (WT = 1)" )
    
    if (!is.null(D.err)){
      D.err <- rbind(D.err, c(NA, 0,  0, 0, NA)) # TODO find a better may to pass sd/error
      D.err = D.err[SORTED,]
      
      pd <- position_dodge(0.1) # move them .05 to the left and right
      
      p <-  p + geom_errorbar(mapping=aes(x=factor(rank( D[,metric], ties.method = "first" )),
                                          ymin=D[,metric] - (D.err[,metric] * ERR.factor), 
                                          ymax=D[,metric] + (D.err[,metric] * ERR.factor), 
                                          colour=factor(Class)), 
                              width=DEFAULT_LN_WIDTH, 
                              size=DEFAULT_LN_SIZE, 
                              position=pd) 
    }
    
    if (length(LETHAL) > 0){
      p  = p + geom_hline(yintercept = LETHAL_PENALTY-0.05, linetype="solid") +
               geom_text(data=data.frame(1),
                  x=80,
                  y=LETHAL_PENALTY,
                  vjust = 1,
                  label = "Lethal", 
                  family = "Arial",
                  fontface = "bold",
                  colour="#696969",
                  #colour="#b32400",
                  size=5)
    }
    
    # gggsave(paste0("wormPointPlot-",metric), DIR=DIR)
    p
    return(p)
  }
# flyPointPlot(D=fly.normalized.df, D.err=fly.err.normalized.df, metric = "Fly.activity", LETHAL = c("4A","E256K"))


# 
# 
# flyAAPlot <- function(D, Y,  DIR="") {
#   # D - Drosophila dataframe
#   
#   D <- merge(D, ClassInformation)  
#   common.genes = D$PTENVariant[!is.na(D$fly.value)]
#   D.common <- D[D$PTENVariant %in% common.genes, ]
#   rownames(D.common) <- D.common$PTENVariant
#   
#   D.common <- D.common[common.genes,]
#   
#   D.AAORDER <- aaOrder(D.common$PTENVariant)
#   
#   D.common <- cbind(D.common, AAORDER=D.AAORDER)
#   
#   g <- guide_legend("title")
#   p <- ggplot(D.common, 
#               aes(x=D.AAORDER, 
#                   y=fly.value, 
#                   fill=factor(D.common$Class)  ))  +
#     
#     geom_point( size = DEFAULT_PT_SIZE*2, shape=21, stat = "identity", color='black', stroke=2) +
#     
#     geom_smooth(data=D.common, size=REGRESSION_LN_SIZE, formula=y~x, method=loess, color='green',aes(x = AAORDER, y = D$fly.value, group=1), se=F) +
#     theme_classic() +
#     geom_hline(aes(yintercept=0), linetype='dotted') +
#     geom_hline(aes(yintercept=1), linetype='dotted') +
#     TITLE_SIZES + 
#     colScaleFill + 
#     labs(title=paste("Amino acid position per estimated variant activity."), 
#          x= "Amino acid position", 
#          y = "Estimated activity (WT = 1)" )  + 
#     scale_fill_manual(name="Class",
#                       values=classColors,
#                       guide = guide_legend(
#                         override.aes = list(
#                           shape = 15,
#                           colour = unique(
#                             as.character(
#                               classColors[unique(sort(D.common$Class))]  ))
#                         ) 
#                       )) 
#   
#   p
#   gggsave("flyYeastVersusAA", DIR=DIR)
#   return(p)
# }


flyYeastAAPlot <- function(D, Y,  DIR="") {
  # D - Drosophila dataframe
  # Y - Yeast dataframe
  # Example: D = fly; Y = yeast.vac14;
  
  common.genes <- intersect(D$PTENVariant[!is.na(D$Fly.activity)], Y$PTENVariant[!is.na(Y$Yeast.activity)])
  
  D.common <- D[D$PTENVariant %in% common.genes, ]
  Y.common <- Y[Y$PTENVariant %in% common.genes, ]
  rownames(Y.common) <- Y.common$PTENVariant 
  rownames(D.common) <- D.common$PTENVariant
  
  Y.common <- Y.common[common.genes,]
  D.common <- D.common[common.genes,]
  
  D.AAORDER <- aaOrder(D.common$PTENVariant)
  Y.AAORDER <- aaOrder(Y.common$PTENVariant)
  
  D.common <- cbind(D.common, AAORDER=D.AAORDER)
  Y.common <- cbind(Y.common, AAORDER=Y.AAORDER)
  
  E.min = sapply(1:length(D.common$Fly.activity),
                 FUN = function(x) return(min(Y.common$Yeast.activity[x],D.common$Fly.activity[x]))
  )
  E.max = sapply(1:length(D.common$Fly.activity),
                 FUN = function(x) return(max(Y.common$Yeast.activity[x],D.common$Fly.activity[x]))
  )
  
  E.colour = as.character(E.min)
  E.colour[T] <- "Yeast > Fly"
  E.colour[Y.common$Yeast.activity < D.common$Fly.activity] <- "Fly > Yeast"
  E.x = D.AAORDER - 1.5
  E.x[Y.common$Yeast.activity < D.common$Fly.activity] = E.x[Y.common$Yeast.activity < D.common$Fly.activity] + (1.5 * 2)
  
  g <- guide_legend("title")
  p <- ggplot(D.common, 
              aes(x=D.AAORDER, 
                  y=Fly.activity, 
                  fill=factor(Y.common$Class)  ))  +
    geom_errorbar(mapping=aes(x=E.x,
                              ymin=E.min,
                              ymax=E.max,
                              colour=E.colour),
                  size = 3.5) +
    scale_colour_manual(name="Comparison", 
                        values=c("green","orange"),
                        guide = guide_legend(
                          override.aes = list( shape = c(15,21) ))
                        ) +
    geom_errorbar(mapping=aes(x=D.AAORDER,
                              ymin=E.min,
                              ymax=E.max),
                  colour='white',
                  size = 2.0) +
    geom_point( size = DEFAULT_PT_SIZE*2, shape=21, stat = "identity", color='black', stroke=2) +
    geom_point( mapping=aes(x=Y.AAORDER, 
                            y=Y.common$Yeast.activity),
                size = DEFAULT_PT_SIZE*2, 
                shape=24, 
                stat = "identity", 
                color='black', 
                stroke=2) +
    theme_classic() + 
    geom_hline(aes(yintercept=0), linetype='dotted') +
    geom_hline(aes(yintercept=1), linetype='dotted') +
    TITLE_SIZES + 
    colScaleFill + 
    labs(title=paste("Amino acid position per estimated variant activity."), 
         x= "Amino acid position", 
         y = "Estimated activity (WT = 1)" )  + 
    scale_fill_manual(name="Class",
                      values=classColors,
                      guide = guide_legend(
                        override.aes = list(
                          shape = c(15,15,15,15,15,15,15,15),
                          colour = unique(
                            as.character(
                              classColors[unique(sort(Y.common$Class))]  ))
                        ) 
                      )) 
  
  p
  
  gggsave("flyYeastVersusAA", DIR=DIR)
  return(p)
}
