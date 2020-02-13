# description: "Project wide constants"
# more: "Use for plotting or logical constants; in any case different set of constants should be swappable. E.g. a constant file could be made to generate plots with finer type, etc."
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
from("utils.R")

## Directory 
PLOT.DIR="plots/"
DATA.DIR="res/"

## Settings
ENABLE_GGGSAVE=FALSE

# Files
PTEN_METADATA = "ExternalData/PTEN.ClassInformation.tsv"
SENTINEL = ""

## Labelling
# Constants
STR_ASD <- "ASD"
STR_EMPTYVECTOR <- "Empty vector"
STR_CLINVAR <- "ClinVar"
STR_WT <- "Wildtype overexpression"
STR_EV <- "Empty vector"
STR_POP <- "ExAC"
STR_HIGHPRED <- "Bioinformatics High impact"
STR_LOWPRED <- "Bioinformatics Low impact"
STR_BIOCHEM <- "Biochemical mutants"

## Rendering
DEFAULT_PT_SIZE <- 3
DEFAULT_LN_WIDTH <- 0
DEFAULT_LN_SIZE <- 1
REGRESSION_LN_SIZE <- 2
DEFAULT_AXIS_TEXT <- 16
DEFAULT_TITLE_TEXT <- 20
DEFAULT_TICK_TEXT <- 14
TITLE_SIZES <- theme(axis.title.x = element_text(size=DEFAULT_AXIS_TEXT),
                     axis.title.y = element_text(size=DEFAULT_AXIS_TEXT),
                     title = element_text(size=DEFAULT_TITLE_TEXT),
                     legend.title = element_text(size=DEFAULT_AXIS_TEXT))
SMALL_TITLE_SIZES <- theme(axis.title.x = element_text(size=DEFAULT_AXIS_TEXT/2),
                     axis.title.y = element_text(size=DEFAULT_AXIS_TEXT/2),
                     title = element_text(size=DEFAULT_TITLE_TEXT/2),
                     legend.title = element_text(size=DEFAULT_AXIS_TEXT/2))

SCATTER_PT_SIZE <- 4
SCATTER_TITLE_SIZES <- TITLE_SIZES + 
  theme(axis.text = element_text(size=DEFAULT_TICK_TEXT))

