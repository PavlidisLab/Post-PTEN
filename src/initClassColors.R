from("src/plotting.R")

# Constants
STR_ASD <- "ASD"
STR_EMPTYVECTOR <- "Control"
STR_CLINVAR <- "ClinVar"
STR_WT <- "Wildtype overexpression"
STR_POP <- "ExAC"
STR_HIGHPRED <- "Bioinformatics High impact"
STR_LOWPRED <- "Bioinformatics Low impact"
STR_BIOCHEM <- "Biochemical mutants"

# Classes
lowPredictedEffectClasses <- c(STR_WT, STR_POP, STR_LOWPRED)
highPredictedEffectClasses <- c(STR_EMPTYVECTOR, STR_ASD, STR_CLINVAR, STR_HIGHPRED)

# Color scheme
classColors <- c(   
  
  "Reduction of function"="#000000", # Teal
  "Rescue"="#808000", # Olive
  "daf-18 + PTEN"="#000000",
  
  "ASD" = "#4971FA", # Blue
  "ASD/ID" = "#4971FA", # Blue    
  
  "Human WT"="#808000", # Olive, equivalent to Rescue
  
  "Biochemical mutants" = "#01FCC4", # Windex blue
  
  "Bioinformatics High impact"="#FF0000", # Red
  "Bioinformatics Low impact"="#FF8FDD",  # Pink
  
  "ClinVar"= "#FF8C00", # Orange
  #"#FEFB06", # Yellow
  # "Cancer"="#FF8C00", # Orange
  
  "ExAC"= "#A020F0",  # Purple ; #FFFFFF/White doesn't work without the border.
  "Population"= "#A020F0",  # Purple ; #FFFFFF/White doesn't work without the border.
  
  "Wildtype overexpression"= "#000000", # Black
  "PTEN overexpression"= "#000000", # Black
  "Wildtype"= "#000000", # White
  
  "Empty vector"="#00FF00", # GREY
  "Control"="black",
  
  "Unclassified"="#B4B4B4",
  "ASD + Cancer"="#a52a2a", # Brown
  
  "Hek.stability"="green",
  "Fly.activity"="orange",
  "Yeast.stability"="pink",
  "Yeast.activity"="red",
  "Total.Dendrite.Length"="purple",
  "Axon.activity"="brown",
  "Worm.activity"="black",
  "PSD.95.Density"="blue",
  
  # "PHTS" = "orange",
  "DD" = "green",
  "Other" = "#a52a2a",
  
  "Variant of Clinical Interest" = "#a52a2a",
  "Other Variants of Clinical Interest" = "#a52a2a",
  "Somatic cancer"="#FF8C00",
  
  "ASD" = "#4971FA", # Blue,
  "Biochemical Control"= "#01FCC4", # Windex blue
  "Cancer"="#FF8C00", # Orange
  "PHTS"="#a52a2a", # Brown
  "Population Control"= "#A020F0",  # Purple  
  "Pred High Impact"="#FF0000", # Red
  "Pred Low Impact"="#FF8FDD"  # Pink
  
)
HIGH_CONTRAST = c("#FEFB06", "#01FCC4", "#FF8FDD") # Point that are high contrast, so white text wouldn't show up.

# Useful for barplots
colScale <-  scale_colour_manual(name="Class", values = classColors)
colScaleFill <-  scale_fill_manual(name="Class", values = classColors)

