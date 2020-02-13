# description: "Plotting rendering, saving, display functions"
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
from("utils.R")

source("src/initClassColors.R")
source("src/external.multiplot.R")

colScale <-  scale_colour_manual(name="Class", values = classColors)
my_palette <- colorRampPalette(c("#000000" ,"#800000" ,"#FF8000" ,"#FFFF00", "#FFFFFF"))(20)
col_breaks = c(seq(0 , 0.33, length=7),  # (optional) defines the color breaks manually for a "skewed" color transition
               seq(0.34 , 0.66,  length=7),
               seq( 0.67 , 1.0  ,  length=7))

CI_FACTOR = 1.96

SP.THEME =  list(
  theme(legend.position="none"),
  labs(title=NULL)
)
SP_PT_SIZE = 2
SP.TITLE_SIZES = SMALL_TITLE_SIZES
SP.PREFIX = "Pearson/Spearman:"


gggsave <- 
  # Save ggplot as postscript/png format to DIR/Prefix.{png,ps}  
  function(Prefix, DIR=PLOT.DIR, scale=1, width=NULL, height=NULL){
    print(paste0( "Saving files at ", DIR, Prefix,".png" ))
    if (ENABLE_GGGSAVE){
      ggsave(paste0( DIR, Prefix,".png" ), scale=scale, width = width, height = height )
      ggsave(paste0( DIR, Prefix,".ps"  ), scale=scale, width = width, height = height )
    }
  }

###################################################################################################################################################

sfari_heatmap <- 
  # A heatmap displaying all the results from the collaboration.
  function(DATAFRAME, 
           Colours=NULL, 
           SortOrder = TRUE, 
           Agg.Func=c("mean"), 
           Agg.Colours=c("black"), 
           Title.Main="Variant overlap", ...) {
    
    MATRIX <- t(
      apply(X = DATAFRAME,
            MARGIN = 2, 
            FUN = function(x) {
              y <- (x - min(x, na.rm = T))
              z <- y / max(y, na.rm = T)
              return(z)
            })
    )
    
    if (is.null(Colours)){
      Colours <- colorRampPalette((c("#000000" ,"#800000" ,"#FF8000" ,"#FFFF00", "#FFFFFF")))(20)
    }
    
    
    
    hm <- main_heatmap(MATRIX,
                       name = "Main<br>Title", 
                       x_categorical = TRUE,
                       y_categorical = TRUE,
                       layout = list(font = list(size = 14)),
                       colors = Colours,
                       ...
    ) %>%
      add_col_labels(ticktext = colnames(MATRIX),
                     font = list(size = 12)) %>%
      add_row_labels(size = 0.3,font = list(size = 12)) %>%
      add_col_title(Title.Main, side= "top") 
    
    
    return(hm)
  }



makeColScale <- function(df, Variant="Variant", Class="Class", type='colour') {
  # Make a color scale from Variants based on the Class in a dataframe
  
  classesVariants<-unique(df[,c(Variant, Class)])
  classes<-classesVariants[,-1,drop=F]
  
  classColorsOrder <- as.numeric(c())
  for (name in names(classColors)) {
    #print(which( (name == classes$Class) == TRUE ))
    rows <- which( (name == classes$Class) == TRUE )
    classColorsOrder <- c( classColorsOrder, rows)
  }
  
  colScale <-  scale_colour_manual(name="Class", values = classColors)
  if (type == "fill"){
    colScale <-  scale_fill_manual(name="Class", values = classColors)
  }
  
  return(colScale)
}
#
PALETTE.correlation=c("#de2d26",
                      "#fc9272",
                      "#FFFFFF",
                      "#a1d99b",
                      "#31a354")
# 
# PALETTE=c("#de2d26",
#           "#fc9272",
#           "#000000",
#           "#74c476",
#           "#31a354",
#           "#006d2c")
# PALETTE = c("#FF0000", "#000000" , "#00FF00")
# PALETTE = c("#FF0000", "#000000" , "#74c476", "#00FF00")
PALETTE.default= c("#ca0020",
                   "#f4a582",
                   "#FFFFFF",
                   "#92c5de",
                   "#0571b0")
#PALETTE = c("#FF0000", "#000000", "#00FF00")


heatmap.common <- 
  # Common data heatmap
  function(DATA,COLOURS=NULL){
    heatmap.2(DATA, 
              main = "Cross-assay comparison", # heat map title
              notecol="black",      # change font color of cell labels to black
              density.info="none",  # turns off density plot inside color legend
              trace="none",         # turns off trace lines inside the heat map
              margins = c(12,30),     # widens margins around plot
              col=colorRampPalette(PALETTE)(20), # use on color palette defined earlier
              #col=my_palette,       # use on color palette defined earlier
              #breaks=col_breaks,    # enable color transition at specified limits
              breaks=seq(from = -0.5, to = 2.5, length.out = 21),    # enable color transition at specified limits
              # dendrogram="row",     # only draw a row dendrogram
              dendrogram='none',
              Colv = F, #ORDER,
              Rowv = F, 
              cexCol = 2.5,
              cexRow = 3,
              na.color = 'darkgrey', 
              labCol = colnames(DATA),
              keysize=1,
              ColSideColors=COLOURS,
              symkey=F
    )   
  }


pheatmap.common <-
  # Common data heatmap
  # Examples:
  #   DATA = AllCommonData.df; COLOURS=MATRIX.CLASSCOLOURS;
  function(DATA,COLOURS=NULL, PALETTE=PALETTE.default, breaks=NULL, 
           CLUSTER_COLS=F,
           CLUSTER_ROWS=F, ...){
    
    annotations = NULL
    anno_colors = NULL
    if (!is.null(COLOURS)){
      annotations <- data.frame(Class = names(COLOURS))
      row.names(annotations) <- colnames(DATA)
      
      Class        <- unique(as.character(COLOURS))
      names(Class) <- unique(names(COLOURS))
      anno_colors <- list(Class = Class)
    }
    
    if (is.null(breaks)){
      # breaks = seq(from = -0.5, to = 2.5, length.out = 21)  
      breaks = seq(from = -0.2, to = 1.2, length.out = 21)  
    }
    
    DATA[ DATA < min(breaks) ] <- min(breaks)
    DATA[ DATA > max(breaks) ] <- max(breaks)
    
    pheatmap.annotations = pheatmapAnnotations(DATA)
    
    
    pheatmap::pheatmap(mat=DATA,
                       labels_row = as.label(row.names(DATA)),
                       color = colorRampPalette(PALETTE)(20),
                       breaks = breaks,    # enable color transition at specified limits
                       cluster_cols = CLUSTER_COLS,
                       cluster_rows = CLUSTER_ROWS,
                       fontsize = 20, 
                       annotation = pheatmap.annotations$annotations,
                       annotation_colors = pheatmap.annotations$colors,
                       #gaps_row = 2,
                       lwid=c(0.1,0.1),
                       lhei=c(0.1,0.1), 
                       na_col = 'black'
                       ,...
    )
  }
# pheatmap.common(AllCommonData.df, COLOURS=MATRIX.CLASSCOLOURS)

pheatmapAnnotations <-
  # Doc
  function(DATA,
           ClassInformation,
           display.class=T,
           display.exac=T,
           display.clinvar=T,
           Index="Genotype"){
    
    # Make colours for variant classes.
    MATRIX.CLASSCOLOURS <- c()
    considered.classes = as.character(unique(ClassInformation$Class))
    
    for (variant in colnames(DATA)){
      current.class = as.character(ClassInformation[variant == ClassInformation[,Index],"Class"])
      if (is_empty(current.class))
        MATRIX.CLASSCOLOURS <- c(MATRIX.CLASSCOLOURS, "#666666") # Unclassified
      else if (variant == "WT") 
        MATRIX.CLASSCOLOURS <- c(MATRIX.CLASSCOLOURS, "#FFFFFF") # WT
      else if (current.class %in% considered.classes ){
        current.colour = classColors[ current.class ]
        MATRIX.CLASSCOLOURS <- c(MATRIX.CLASSCOLOURS, current.colour)
      }
    }
    
    names(MATRIX.CLASSCOLOURS)[names(MATRIX.CLASSCOLOURS) == ""] = "Unclassified"
    if (!is.null(MATRIX.CLASSCOLOURS)){
      annotations <- data.frame(Class = names(MATRIX.CLASSCOLOURS))
      row.names(annotations) <- colnames(DATA)
      
      Class        <- unique(as.character(MATRIX.CLASSCOLOURS))
      names(Class) <- unique(names(MATRIX.CLASSCOLOURS))
      anno_colors <- list(Class = Class)
    }
    
    CLINVAR = annotations.clinvar[ , c("PTENVariant",  "CLNSIG") ];
    if (!is.null(CLINVAR)){
      CLINVAR$CLNSIG <- gsub(CLINVAR$CLNSIG, pattern = "/", replacement = ".")
      # CLINVAR$CLNSIG <-
      sapply(1:nrow(CLINVAR), FUN = function(I){
        VARIANT = CLINVAR$PTENVariant[I]
        noquote( paste(sep = "|", (unique((CLINVAR[CLINVAR$PTENVariant == VARIANT, "CLNSIG"]))) ) )
      })
      CLINVAR.dedup = CLINVAR[ !duplicated(CLINVAR$PTENVariant), ]
      CLINVAR.CLASSES = CLINVAR.dedup[ match(rownames(annotations), CLINVAR.dedup$PTENVariant), ]
      CLINVAR.CLASSES[is.na(CLINVAR.CLASSES$CLNSIG), "CLNSIG"] <- "NA"
      
      # ClinVar multiple reports
      
      CLINVAR.CLASSES$CLNSIG[
        grepl(CLINVAR.CLASSES$CLNSIG, pattern = "Conflicting", ignore.case = T)] <- 'Conflicting'# Important to do this one before "Conflicting reports of *pathogenic*ity"
      
      CLINVAR.CLASSES$CLNSIG[
        grepl(CLINVAR.CLASSES$CLNSIG, pattern = "Pathogenic", ignore.case = T) &
          grepl(CLINVAR.CLASSES$CLNSIG, pattern = "Benign", ignore.case = T)] <- 'Conflicting'
      
      CLINVAR.CLASSES$CLNSIG[
        grepl(CLINVAR.CLASSES$CLNSIG, pattern = "Uncertain", ignore.case = T) |
          grepl(CLINVAR.CLASSES$CLNSIG, pattern = "Unknown", ignore.case = T)] <- 'Unknown'
      
      CLINVAR.CLASSES$CLNSIG[
        grepl(CLINVAR.CLASSES$CLNSIG, pattern = "Benign", ignore.case = T)] <- 'Benign'
      
      CLINVAR.CLASSES$CLNSIG[
        grepl(CLINVAR.CLASSES$CLNSIG, pattern = "Pathogenic", ignore.case = T)] <- 'Pathogenic'
      
      
      annotations$ClinVar = CLINVAR.CLASSES$CLNSIG
      
      
      CLNSIG        <- colorRampPalette(c("red", "white"))( nrow(unique(CLINVAR.CLASSES[, "CLNSIG"])) )
      names(CLNSIG) <- as.character(t(unique(CLINVAR.CLASSES[, "CLNSIG"])))
      
      CLNSIG["NA"] <- 'black'
      
      CLNSIG[grepl(names(CLNSIG), pattern = "Benign", ignore.case = T) ] <- '#0000FF'
      CLNSIG[grepl(names(CLNSIG), pattern = "Pathogenic", ignore.case = T) ] <- '#ff6666' # Important to do this one before "Conflicting reports of *pathogenic*ity"
      
      CLNSIG[grepl(names(CLNSIG), pattern = "Uncertain", ignore.case = T) ] <- '#D3D3D3'
      CLNSIG[grepl(names(CLNSIG), pattern = "Conflicting", ignore.case = T) ] <- 'cadetblue'
      CLNSIG[grepl(names(CLNSIG), pattern = "Unknown", ignore.case = T) ] <- '#D3D3D3'
      
      anno_colors$`ClinVar` = CLNSIG
    }
    
    POPULATION.FREQ = annotations.annovar[, c("PTENVariant", "ExAC_ALL")] # TODO: This was just pasted here to make it work, but really it should be set elsewhere.
    if (!is.null(POPULATION.FREQ)){
      POPULATION.FREQ$ExAC_ALL <- sapply(1:nrow(POPULATION.FREQ),
                                         FUN = function(I){
                                           VARIANT = POPULATION.FREQ$PTENVariant[I]
                                           sum(unique(POPULATION.FREQ[POPULATION.FREQ$PTENVariant == VARIANT, "ExAC_ALL"]), na.rm = T)
                                         })
      POPULATION.FREQ.dedup = POPULATION.FREQ[ !duplicated(POPULATION.FREQ$PTENVariant), ]
      POPULATION.FREQ.CLASSES = POPULATION.FREQ.dedup[ match(rownames(annotations), POPULATION.FREQ.dedup$PTENVariant), ]
      POPULATION.FREQ.CLASSES$ExAC_ALL[is.na(POPULATION.FREQ.CLASSES$ExAC_ALL) | POPULATION.FREQ.CLASSES$ExAC_ALL == 0] <- "N = 0"
      POPULATION.FREQ.CLASSES$ExAC_ALL[POPULATION.FREQ.CLASSES$ExAC_ALL != "N = 0" ] <- "N >= 1"
      
      annotations$`ExAC` = POPULATION.FREQ.CLASSES$ExAC_ALL
      
      
      ExAC_ALL        <- colorRampPalette(c("red", "white"))( nrow(unique(POPULATION.FREQ.CLASSES[, "ExAC_ALL"])) )
      names(ExAC_ALL) <- as.character(t(unique(POPULATION.FREQ.CLASSES[, "ExAC_ALL"])))
      
      ExAC_ALL[ "N >= 1" ] <- '#cc99ff'
      ExAC_ALL["N = 0"] <- 'black'
      
      anno_colors$`ExAC` = ExAC_ALL
    }
    annotations$PHTS = if_else(ASDFirstClass[match(colnames(DATA), ASDFirstClass$PTENVariant), "PHTS.variant"] == "1",
                               "Yes",
                               "No")
    anno_colors$PHTS = c(Yes="#cc6666", No='black')
    
    annotations$`Somatic cancer` = if_else(ASDFirstClass[match(colnames(DATA), ASDFirstClass$PTENVariant), "Cancer.variant"] == "1",
                                           "Yes",
                                           "No")
    anno_colors$`Somatic cancer` = c(Yes="#ffcc99", No='black')
    
    return(list(annotations=annotations, colors=anno_colors))
  }

pheatmapAnnotations.MEM <-
  # Doc
  function(DATA,
           ClassInformation,
           display.class=T,
           display.exac=T,
           display.clinvar=T,
           Index="Genotype"){
    
    # Make colours for variant classes.
    MATRIX.CLASSCOLOURS <- c()
    considered.classes = as.character(unique(ClassInformation$Class))
    
    # CLASS COLOURS
    for (variant in colnames(DATA)){
      current.class = as.character(ClassInformation[variant == ClassInformation[,Index],"Class"])
      if (is_empty(current.class))
        MATRIX.CLASSCOLOURS <- c(MATRIX.CLASSCOLOURS, "#666666") # Unclassified
      else if (variant == "WT") 
        MATRIX.CLASSCOLOURS <- c(MATRIX.CLASSCOLOURS, "#FFFFFF") # WT
      else if (current.class %in% considered.classes ){
        current.colour = classColors[ current.class ]
        MATRIX.CLASSCOLOURS <- c(MATRIX.CLASSCOLOURS, current.colour)
      }
    }
    
    names(MATRIX.CLASSCOLOURS)[names(MATRIX.CLASSCOLOURS) == ""] = "Unclassified"
    if (!is.null(MATRIX.CLASSCOLOURS)){
      annotations <- data.frame(Class = names(MATRIX.CLASSCOLOURS))
      row.names(annotations) <- colnames(DATA)
      
      Class        <- unique(as.character(MATRIX.CLASSCOLOURS))
      names(Class) <- unique(names(MATRIX.CLASSCOLOURS))
      anno_colors <- list(Class = Class)
    }
    
    # DENOVO STATUS (De.Novo)
    ## Classes
    annotations$De.Novo <- as.character(ClassInformation$De.Novo[ match(rownames(annotations), ClassInformation$Genotype) ])
    annotations$De.Novo[is.na(annotations$De.Novo) | annotations$De.Novo == ""] <- "N/A"
    annotations$De.Novo
    ## Colours
    denovo.colors = c("#ff60fc", "black", "white")
    names(denovo.colors) <- c("Y", "N", "N/A")
    anno_colors$De.Novo <- denovo.colors
    
    # ASD STATUS (ASD.Associated)
    ## Classes
    annotations$ASD.Associated <- as.character(ClassInformation$ASD.Associated[ match(rownames(annotations), ClassInformation$Genotype) ])
    annotations$ASD.Associated[is.na(annotations$ASD.Associated) | annotations$ASD.Associated == ""] <- "N/A"
    annotations$ASD.Associated
    ## Colours
    ASD.colors = c("#00177f", "black", "white")
    names(ASD.colors) <- c("Y", "N", "N/A")
    anno_colors$ASD.Associated <- ASD.colors
    
    #DD STATUS (Developmental.Delay)
    ## Classes
    annotations$DD.Associated <- as.character(ClassInformation$Developmental.Delay[ match(rownames(annotations), ClassInformation$Genotype) ])
    annotations$DD.Associated[is.na(annotations$DD.Associated) | annotations$DD.Associated == ""] <- "N/A"
    annotations$DD.Associated
    ## Colours
    DD.colors = c("#007f17", "black", "white")
    names(DD.colors) <- c("Y", "N", "N/A")
    anno_colors$DD.Associated <- DD.colors

    # CANCER STATUS (Somatic.Cancer)
    ## Classes
    annotations$Somatic.Cancer <- as.character(ClassInformation$Somatic.Cancer[ match(rownames(annotations), ClassInformation$Genotype) ])
    annotations$Somatic.Cancer[is.na(annotations$Somatic.Cancer) | annotations$Somatic.Cancer == ""] <- "N/A"
    annotations$Somatic.Cancer
    ## Colours
    Somatic.Cancer.colors = c("#cc6a10", "black", "white")
    names(Somatic.Cancer.colors) <- c("Y", "N", "N/A")
    anno_colors$Somatic.Cancer <- Somatic.Cancer.colors
    
    # PHTS STATUS (PHTS)
    ## Classes
    annotations$PHTS <- as.character(ClassInformation$PHTS[ match(rownames(annotations), ClassInformation$Genotype) ])
    annotations$PHTS[is.na(annotations$PHTS) | annotations$PHTS == ""] <- "N/A"
    annotations$PHTS
    ## Colours
    PHTS.colors = c("#774a20", "black", "white")
    names(PHTS.colors) <- c("Y", "N", "N/A")
    anno_colors$PHTS <- PHTS.colors
    
    annotations = annotations[, c(setdiff(names(annotations), "Class"), "Class")] # Force Class to be at the top.
    return(list(annotations=annotations, colors=anno_colors))
  }



pheatmap.long.MEM <-
  # Long form heatmap
  function(DATA, 
           ClassInformation,
           COLOURS=NULL, 
           PALETTE=PALETTE.default, 
           breaks=NULL,
           FONT_SIZE=8,
           X_FONT_SIZE=NULL,
           Y_FONT_SIZE=NULL,
           SORTBY="Class",
           gaps_col=NULL,
           METRICS=NULL,
           METRICS.subsort=NULL,
           showWT=F,
           showMissenseOnly = T,
           showTruncatingOnly = F, 
           CLINVAR=NULL,
           POPULATION.FREQ=NULL,
           Index="Genotype",
           WT.VARIANT = "WT",
           cellheight = 8,
           cellwidth = 7,
           ...){
    
    if (is.null(X_FONT_SIZE)){ X_FONT_SIZE = FONT_SIZE }
    if (is.null(Y_FONT_SIZE)){ Y_FONT_SIZE = FONT_SIZE }
    
    if(showMissenseOnly){
      DATA = missenseOnly(DATA, Index = Index)
    } else if (showTruncatingOnly){
      DATA = truncatingOnly(DATA, Index = Index)
    } else if (showMissenseOnly == T & showTruncatingOnly == T){
      print("[ERROR] Cannot be missense only and truncating only.")
      exit(-1)
    }
    
    if (!showWT){
      DATA = DATA[DATA[,Index] != WT.VARIANT,]
    }
    
    annotations = NULL
    anno_colors = NULL
    
    if (is.null(breaks)){
      # Default heatmap value range ([-0.2 , 1.2])
      breaks = seq(from = -0.2, to = 1.2, length.out = 21)  
    }
    
    
    DATA[ DATA < min(breaks) ] <- min(breaks)
    DATA[ DATA > max(breaks) ] <- max(breaks)
    
    if (is.null(METRICS)){
      # TODO: Update scores
      METRICS=c("AssayMean",
                "Hek.activity",
                "Yeast.activity",
                "Fly.activity",
                "mighell_bin",
                "Yeast.stability",
                "Hek.stability",
                "matreyek_score",
                "SNAP2",
                "CADD13Raw")
    }
    rownames(DATA) <- as.character(DATA[,Index])
    
    AVERAGE_SCORE_ORDER = rank(rowMeans(DATA[,METRICS.subsort], na.rm = T), ties.method = 'first')
    
    if ( SORTBY == "mean" ) {
      # Sort by mean score
      SORT_ORDER = AVERAGE_SCORE_ORDER
      DATA = DATA[SORT_ORDER, METRICS]
      
    } else if (SORTBY == "Class") {
      DATA$ClassSorting <- factor(DATA$Class, 
                                  c("Biochemical Control",
                                    "ASD",
                                    "PHTS",
                                    "Cancer",
                                    "Pred High Impact",
                                    "Pred Low Impact",
                                    "Population Control"))
      
      
      SORT_ORDER = with(DATA, order(DATA$ClassSorting, AVERAGE_SCORE_ORDER) )
      DATA = DATA[SORT_ORDER, ]
      
      if ( is.null(gaps_col) ) {
        # Set column gaps based on Class Sorting
        gaps_col = sapply( 
          1:length(unique((DATA$ClassSorting))),
          FUN = function(x)
            sum(as.numeric(table(DATA$ClassSorting)[unique(DATA$ClassSorting)][1:x]))
        )
      }
      DATA = DATA[, METRICS]
    }  else {
      # Sort by column specified in SORTBY
      DATA = DATA[order(DATA[,SORTBY]),METRICS]
    }
    
    DATA = t(DATA[,setdiff(names(DATA), c(Index, "Class"))])
    pheatmap.annotations = pheatmapAnnotations.MEM(DATA, ClassInformation)
    pheatmap::pheatmap(mat=DATA,
                       labels_row = as.label(row.names(DATA)),
                       color = colorRampPalette(PALETTE)(20),
                       breaks = breaks,    # enable color transition at specified limits
                       cluster_cols = F,
                       cluster_rows = F,
                       fontsize_col = X_FONT_SIZE,
                       fontsize_row = Y_FONT_SIZE,
                       annotation = pheatmap.annotations$annotations,
                       annotation_colors = pheatmap.annotations$colors,
                       #gaps_row = 2,
                       na_col = 'black',
                       gaps_col = gaps_col[!is.na(gaps_col)], 
                       lwid=c(0.1,0.1),
                       lhei=c(0.1,0.1),
                       cellheight = cellheight,
                       cellwidth = cellwidth,
                       # gaps_row=gaps_row#,
                       ...
                       
    )
  }

pheatmap.long <-
  # Long form heatmap
  function(DATA, 
           COLOURS=NULL, 
           PALETTE=PALETTE.default, 
           breaks=NULL,
           FONT_SIZE=8,
           X_FONT_SIZE=NULL,
           Y_FONT_SIZE=NULL,
           SORTBY="Class",
           gaps_col=NULL,
           METRICS=NULL,
           showWT=F,
           showMissenseOnly = T,
           showTruncatingOnly = F, 
           CLINVAR=NULL,
           POPULATION.FREQ=NULL,  ...){
    
    if (is.null(X_FONT_SIZE)){ X_FONT_SIZE = FONT_SIZE }
    if (is.null(Y_FONT_SIZE)){ Y_FONT_SIZE = FONT_SIZE }
    
    if(showMissenseOnly){
      DATA = missenseOnly(DATA)
    } else if (showTruncatingOnly){
      DATA = truncatingOnly(DATA)
    } else if (showMissenseOnly == T & showTruncatingOnly == T){
      print("[ERROR] Cannot be missense only and truncating only.")
      exit(-1)
    }
    
    if (!showWT){
      DATA = DATA[DATA$PTENVariant != "WT",]
    }
    
    annotations = NULL
    anno_colors = NULL
    
    if (is.null(breaks)){
      # Default heatmap value range ([-0.2 , 1.2])
      breaks = seq(from = -0.2, to = 1.2, length.out = 21)  
    }
    
    
    DATA[ DATA < min(breaks) ] <- min(breaks)
    DATA[ DATA > max(breaks) ] <- max(breaks)
    
    if (is.null(METRICS)){
      METRICS=c("AssayMean",
                "Hek.activity",
                "Yeast.activity",
                "Fly.activity",
                "mighell_bin",
                "Yeast.stability",
                "Hek.stability",
                "matreyek_score",
                "SNAP2",
                "CADD13Raw")
    }
    
    rownames(DATA) <- DATA$PTENVariant
    
    if ( length(setdiff(c("Yeast.activity", "Fly.activity"), METRICS)) == 0 ) { # Yeast/Fly are in metrics
      ACTIVITY_ORDER =
        sapply( 1:nrow(DATA),
                FUN = function(x){
                  mean(
                    as.numeric(
                      DATA[ x, c("Yeast.activity", "Fly.activity") ]
                    ), na.rm=T)
                })
    } else {
      # Keep order the same.
      ACTIVITY_ORDER =  sapply( 1:nrow(DATA),
                                FUN = function(x){
                                  mean(
                                    as.numeric(
                                      DATA[ x, ]
                                    ), na.rm=T)
                                })
    }
    
    
    if ( SORTBY == "mean" ) {
      # Sort by mean score
      SORT_ORDER = order(ACTIVITY_ORDER) 
      DATA = DATA[SORT_ORDER, METRICS]
      
    } else if (SORTBY == "ClassASD") {
      DATA$ClassSorting <- factor(DATA$Class, 
                                  c("Biochemical mutants",
                                    "ASD",
                                    # "ASD + Cancer",
                                    "Variant of Clinical Interest",
                                    "Cancer",
                                    "Bioinformatics High impact",
                                    "Bioinformatics Low impact",
                                    "ExAC"))
      
      
      SORT_ORDER = with(DATA, order(DATA$ClassSorting, ACTIVITY_ORDER) )
      DATA = DATA[SORT_ORDER, ]
      
      if ( is.null(gaps_col) ) {
        gaps_col = sapply(1:length(unique((DATA$ClassSorting))), FUN = function(x) sum( as.numeric(table(DATA$ClassSorting)[unique(DATA$ClassSorting)][1:x] )))
        
      }
      DATA = DATA[, METRICS]
      
    }  else if (SORTBY == "Class") {
      DATA$Class <- factor(DATA$Class, 
                           c("Biochemical mutants",
                             "ASD",
                             #"ASD + Cancer",
                             # "Variant of Clinical Interest",
                             "Other Variants of Clinical Interest",
                             # "Cancer",
                             "Somatic cancer",
                             "Bioinformatics High impact",
                             "Bioinformatics Low impact",
                             "ExAC"))
      
      DATA$Class[ DATA$Class == "ASD + Cancer"] <- factor("ASD")
      
      SORT_ORDER = with(DATA, order(DATA$Class, ACTIVITY_ORDER) )
      DATA = DATA[SORT_ORDER, ]
      
      if ( is.null(gaps_col) ) {
        gaps_col = sapply(
          1:length(unique((DATA$Class))), 
          FUN = function(x) 
            sum( as.numeric(table(DATA$Class)[unique(DATA$Class)][1:x] ))
        )
        
      }
      DATA = DATA[, METRICS]
      
    } else {
      DATA = DATA[order(DATA[,SORTBY]),METRICS]
    }
    
    DATA = t(DATA)
    
    pheatmap.annotations = pheatmapAnnotations(DATA)
    
    pheatmap::pheatmap(mat=DATA,
                       labels_row = as.label(row.names(DATA)),
                       color = colorRampPalette(PALETTE)(20),
                       breaks = breaks,    # enable color transition at specified limits
                       cluster_cols = F,
                       cluster_rows = F,
                       fontsize_col = X_FONT_SIZE,
                       fontsize_row = Y_FONT_SIZE,
                       annotation = pheatmap.annotations$annotations,
                       annotation_colors = pheatmap.annotations$colors,
                       #gaps_row = 2,
                       na_col = 'black',
                       gaps_col = gaps_col, 
                       lwid=c(0.1,0.1),
                       lhei=c(0.1,0.1), 
                       ...
                       
    )
  }


MultiOrganismPointPlot <-
  # 
  #
  # Example:
  # M = combined ; M.err = combined.err
  function(M , M.err = NULL, METRICS=NULL, CI=TRUE, DIR=""){
    M$Class <- as.factor(M$Class)
    colScale <- makeColScale(M, Variant="PTENVariant")  
    
    # Set minimum to 0.0
    row.names(M) <- NULL
    
    if (is.null(METRICS)){
      METRICS = c("PTENVariant",
                  "PSD.95.Density",
                  "Total.Dendrite.Length",
                  "Yeast.activity",
                  "Axon.activity",
                  "Fly.activity",
                  "Worm.activity",
                  "Hek.stability",
                  "Yeast.stability",
                  "Class")
    }
    
    VLINE_IDX = min(match(c("Hek.stability", "Yeast.stability"), METRICS)) - 1.5 # Hack; assumes hek is the first protein index is listed here.
    
    
    V = melt(M[,METRICS], id.vars = c("PTENVariant","Class"))
    
    ## Handle (optional) error bars
    ERRORBARS = NULL
    if (!is.null(M.err)) {
      # Handle error bar transformations  
      M.err$Class <- as.factor(M.err$Class)
      row.names(M.err) <- NULL
      METRIC.err = intersect(names(M.err),METRICS)
      
      E = merge(V[,c("PTENVariant", "variable")],
                melt(M.err[,METRIC.err], id.vars = c("PTENVariant","Class")), 
                by=c("PTENVariant", "variable"), all.x=F,all.y=T )
      V = merge(V,
                melt(M.err[,METRIC.err], id.vars = c("PTENVariant","Class"))[,c("PTENVariant", "variable")], 
                by=c("PTENVariant", "variable"), all.x=T,all.y=F )
      
      E.filter = V$variable %in% METRIC.err
      
      # summary(E)
      # summary(V[E.filter,])
      # View(E)
      # View(V[E.filter,])
      
      COMPARISON = compare(E,V[E.filter,] )
      
      stopifnot(
        all(
          COMPARISON[2]$detailedResult[c(1,2,3)]
        )
      )
      # print("Compare passed.")
      
      if (CI) {
        # Set STD.ERR to Confidence intervals
        E$value[ E$variable %in% c("Axon.activity","Fly.activity","PSD.95.Density", "Total.Dendrite.Length", "Worm.activity","Worm.activity","Yeast.activity")] <- 
          E$value[ E$variable %in% c("Axon.activity","Fly.activity","PSD.95.Density", "Total.Dendrite.Length", "Worm.activity","Worm.activity","Yeast.activity")] * CI_FACTOR
      }
      E.upper = V$value[E.filter] + E$value
      E.lower = V$value[E.filter] - E$value
      
      ERRORBARS = geom_errorbar(data=E,
                                mapping = aes(x = variable, ymin = E.lower, ymax = E.upper, colour=Class),
                                stat = "identity",
                                width=DEFAULT_LN_WIDTH, 
                                size=DEFAULT_LN_SIZE)
      
    }
    
    p <- ggplot(V) + 
      ERRORBARS +
      geom_point(data=V,
                 mapping = aes(x = variable, y = value, colour=Class),
                 size = DEFAULT_PT_SIZE*2, 
                 stat = "identity") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) +
      colScale + 
      TITLE_SIZES +
      geom_vline(xintercept = VLINE_IDX, linetype="dotted") + # Protein start
      #geom_vline(xintercept = 8.5, linetype="dotted") + # Protein end
      geom_hline(yintercept = 1, linetype="dashed") +
      geom_hline(yintercept = 0, linetype="dashed") +
      labs(title=paste0("Cross-assay estimated variant activity"), 
           x= "Variants", 
           y = "Estimated activity (WT = 1)" ) +
      facet_wrap( ~ PTENVariant) +
      theme(strip.text.x = element_text(size = 20))
    
    return(p)
  }

theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold",size = rel(1)),
      axis.title.y = element_text(angle=90,vjust =2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2),
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_blank(), #element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size= unit(0.2, "cm"),
      legend.margin = unit(0, "cm"),
      legend.title = element_text(face="italic"),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
    ))
  
}

theme_scatterpairs <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  theme(strip.placement = "outside") +
    (theme_foundation(base_size=14, base_family="helvetica")
     + theme(
       plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
       text = element_text(),
       panel.background = element_rect(colour = NA),
       plot.background = element_rect(colour = NA),
       axis.title = element_text(face = "bold",size = rel(1)),
       axis.title.y = element_text(angle=90,vjust =2),
       axis.title.x = element_text(vjust = -0.2),
       axis.text = element_text(),
       axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2),
       axis.line = element_line(colour="black"),
       axis.ticks = element_line(),
       panel.grid.major = element_blank(), #element_line(colour="#f0f0f0"),
       panel.grid.minor = element_blank(),
       legend.key = element_rect(colour = NA),
       legend.position = "bottom",
       legend.direction = "horizontal",
       legend.key.size= unit(0.5, "cm"),
       legend.margin = unit(0, "cm"),
       legend.title = element_text(face="italic"),
       plot.margin=unit(c(10,5,5,5),"mm"),
       strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
       strip.text = element_text(face="bold")
     )
    )
  
}

get_legend <- 
  # Only return the legend of a ggplot plot.
  ## Sources:
  ## [1] https://stackoverflow.com/questions/11883844/inserting-a-table-under-the-legend-in-a-ggplot2-histogram
  ## [2] https://stackoverflow.com/questions/12041042/how-to-plot-just-the-legends-in-ggplot2
  function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  
  
  grid.newpage()
  grid.draw(legend) 
} 


plotAdjustedActivity <-
  # A quick and dirty visualization of the batch correction effect.
  #
  # means : name of the mean vector
  # sds : name of the sd vector, if not sd.<means>"
  # adjs : name of the sd vector, if not <means>.adj"
  # adjs.sds : name of the sd vector, if not sd.<adjs>"

  function(DF, means, sds=NULL, adjs=NULL, adjs.sds=NULL) {
    if ( is.null(sds) ){
      sds = paste0("sd.", means)
    }
    if ( is.null(adjs) ){
      adjs = paste0(means, ".adj")
    }
    if ( is.null(adjs.sds) ){
      adjs.sds = paste0("sd.", adjs)
    }
    
    print(means)
    print(sds)
    print(adjs)
    print(adjs.sds)
    
    ggplot(DF) + 
      theme_Publication() +
      geom_point(mapping = aes(x = Genotype,
                               y = DF[,means]),
                 color="red") +
      geom_errorbar(mapping = aes(x = Genotype,
                                  ymin = DF[,means] -DF[,sds],
                                  ymax = DF[,means] +DF[,sds]),
                    color="red") + # Plot adjusted
      geom_point(mapping = aes(x = Genotype,
                               y = Chemotaxis_Index.adj),
                 color="green") +
      geom_errorbar(mapping = aes(x = Genotype,
                                  ymin = DF[,adjs] -DF[,adjs.sds],
                                  ymax = DF[,adjs] +DF[,adjs.sds]),
                    color="green")
  }


genericPointPlot <-
  function(DF,
           Value,
           Error,
           Title,
           Subtitle,
           Color.By="Class",
           Group.By=NULL,
           Index="Genotype",
           XLAB="Variants",
           YLAB="Activity"){
    
  DF = DF[!is.na(DF[,Value]) ,]
  
  # Fix colors for NAs
  DF[,Color.By] <- as.character(DF[,Color.By] )
  DF[,Color.By][is.na(DF[,Color.By])] <- "Unclassified"
  DF[DF[,Index] %in% c("WT", "Hum_WT", "cePTEN(rf) + WT PTEN"), Color.By] <- "Control"
  DF[DF[,Index] %in% c("Control", "pEGH","attp2", "GFP", "cePTEN(rf)", "GFP_Control") , Color.By] <- "Control"
  
  if (is.null(Group.By)) {
    DF = DF[order(rank(as.numeric(DF[,Value]), ties.method = "first")),]
    Group.By = Index
  } else if (Group.By == "Class"){
    DF = DF[order(rank(as.numeric(DF[,Value]), ties.method = "first")),]
    DF = DF[order(rank(as.numeric(DF[,Group.By]), ties.method = "first")),]
  } else {
    quit("[ERROR] Unknown Group.By in genericPointPlot")
  }
  DF[,Index] = droplevels(DF[,Index])
  DF[,Index] = factor(DF[,Index], as.character(DF[,Index]))
  
  ggplot(DF) +
      theme_Publication() +
      geom_hline(yintercept = 0.0) +
      geom_hline(yintercept = 1.0) +
      geom_point(mapping = aes(y = DF[,Value],
                               x = DF[,Group.By],
                               color=DF[,Color.By]),
                 size=DEFAULT_PT_SIZE) +
      geom_errorbar(mapping = aes(x=DF[,Group.By],
                                  ymin=DF[,Value]-DF[,Error],
                                  ymax=DF[,Value]+DF[,Error],
                                  color=DF[,Color.By]),
                    width=0,
                    size=1) +
      xlab(XLAB) +
      ylab(YLAB) +
      ggtitle(label = Title, subtitle = Subtitle) +
      theme(axis.text.x = element_text(angle = 90)) + 
      colScale 
}
           

source("src/rat.plotting.R")
source("src/yeast.plotting.R")
source("src/hek.plotting.R")
source("src/worm.plotting.R")
source("src/fly.plotting.R")
source("src/axon.plotting.R")
source("src/comparison.plotting.R")
