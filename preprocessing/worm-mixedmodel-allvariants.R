# description: "Worm data linear model/preprocessing from the raw data."
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "31/01/2019"
source('utils.R')
MODEL = as.formula("Chemotaxis_Index ~ Genotype + (1|Day)") # It's possible here to also use weights=sqrt(Number_of_Animals)
# Note: chemotaxis index = (Number_on_Salt_Side - Number_on_Control_Side)/Number_of_Animals

worm.raw<-read.csv(sep="\t",
                   here::here("RawData/WormChemotaxis/PTEN-Worm.tsv"), 
                   header=1)

worm.raw <- worm.raw[,colSums(!is.na(worm.raw)) > 0] ## Get rid of empty cols
dr<-as.tibble(worm.raw)

PTEN_WT = "cePTEN(rf) + WT PTEN"
PTEN_NULL = "cePTEN(rf)"

# Identify controls
dr$WT<- dr$Genotype == PTEN_WT# FIXME: In the analysis the function uses cePTEN(rf) + WT PTEN to normalize to 1.0
dr$PTENNull <- dr$Genotype == PTEN_NULL
dr$Control <- dr$PTENNull | dr$WT

# Relevel to have genotypes as the background group.
dr$Genotype <- factor(dr$Genotype, c( setdiff(levels(dr$Genotype), c(PTEN_NULL, PTEN_WT) ), c(PTEN_NULL, PTEN_WT) ) )

worms.adjusted.means.df = NULL
dr$Genotype = droplevels(dr$Genotype)
# Apply mixed effect model; treat day as a mixed effect

# Compute p-values from WT and NULL comparisons
dr$Genotype <- factor(dr$Genotype,
                      c( PTEN_NULL, PTEN_WT, setdiff(as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL))) )
modelNULL = lmerTest::lmer(data=dr, MODEL)

dr$Genotype <- factor(dr$Genotype,
                      c( PTEN_WT,PTEN_NULL, setdiff(as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL))) )
modelWT = lmerTest::lmer(data=dr, MODEL)

model = lmerTest::lmer(data=dr, MODEL)
model.rnef = lme4::ranef(model)


pvals.WT = coefficients(summary(modelWT))[,5]
pvals.NULL = coefficients(summary(modelNULL))[,5]
pvals = merge(data.frame(Genotype=names(pvals.WT), pvalue.from.WT=pvals.WT), data.frame(Genotype=names(pvals.NULL), pvalue.from.NULL=pvals.NULL), all.x=T, all.y=T)
###


## Add coefficients and CI
pvals <- merge(pvals, data.frame(Genotype=names(fixef(model)), Coefficients=fixef(model)), all.x=T, all.y=F)
model.confint =  data.frame( model.ci=confint(model) )
colnames(model.confint) <- paste0("Chemotaxis_Index",c("95CI.lo", "95CI.high"))
model.confint$Genotype = rownames(model.confint)

pvals <- merge(pvals, model.confint[grepl(rownames(model.confint), pattern = "Genotype"),], by="Genotype", all.x=T, all.y=T)
pvals$Genotype <- gsub(pvals$Genotype, pattern = "Genotype", replacement = "")

adjusted.means = adjust.data.for.batch.lme4(dr, model, 
                                            model.rnef = model.rnef, 
                                            rnef.names = c("Day"),
                                            # filter.agg = "Chemotaxis_Index",
                                            value.name="Chemotaxis_Index", 
                                            pattern = "Genotype")

adjusted.means <- merge(adjusted.means, pvals, by="Genotype", all.x=T, all.y=F)

adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "Chemotaxis_Index", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)
adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "Chemotaxis_Index.adj", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)
#################################################################################################################################

writeObject(adjusted.means, filename = "worms.MEM.adjusted.means.df", row.names = F)
