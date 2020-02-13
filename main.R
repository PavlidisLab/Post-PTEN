RESET_UTILS=TRUE
source("utils.R")

# Clean-up input data
source("preprocessing/axon-tidy.R")
source("preprocessing/fly-tidy.R")
source("preprocessing/rat-tidy.R")
source("preprocessing/yeast-preprocess-MEM.R")
# Note: Nothing to do for worm preprocessing.

## Estimate activity 
source("preprocessing/axon-mixedmodel-allvariants.R")
source("preprocessing/fly-mixedmodel-allvariants-TDGY.R")
source("preprocessing/rat-mixedmodel-allvariants.R")
source("preprocessing/worm-mixedmodel-allvariants.R")

## Yeast, for each sentinel.
for (YEAST.SENTINEL in c("FIG4", "VAC14", "VAC7", "VAM3", "VAM7", "VPS30", "VPS38", "YPT7")){
  source("preprocessing/yeast-mixedmodel-allvariants.R")
}

