# description: "Axon data mixed model/preprocessing from the raw data using all variants in one model."
# author: "Manuel Belmadani"
# date: "20/02/2019"

source('utils.R')
MODEL = as.formula("axon ~ Genotype + (1 | day) + (1 | plate)")

axon.raw<-read.csv(sep="\t",
                   here::here("Data/PTEN-Axon.tsv"), 
                   header=1)

axon.raw$variant <- factor(axon.raw$variant)
axon.raw$plate = paste0(as.character(axon.raw$variant), "_", as.character(axon.raw$day), "_", as.character(axon.raw$plate))
# axon.raw$plate = paste0(as.character(axon.raw$day), "_", as.character(axon.raw$plate))

axon.raw$day <- factor(axon.raw$day)
axon.raw$plate <- as.factor(axon.raw$plate)

dr <- as.tibble(axon.raw)
dr
names(dr)[1] = "Genotype"

PTEN_WT = "WT"
PTEN_NULL = "GFP"
dr$WT<- dr$Genotype == PTEN_WT
dr$PTEN_NULL <- dr$Genotype == PTEN_NULL
dr$Control <- dr$PTEN_NULL | dr$WT
# dr$axon.log = log10(dr$axon)

# Compute p-values from WT and NULL comparisons
dr$Genotype <- factor(dr$Genotype,
                      c( PTEN_NULL, PTEN_WT, setdiff(as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL))) )
modelGFP = lmerTest::lmer(data=dr, MODEL)
dr$Genotype <- factor(dr$Genotype,
                      c( PTEN_WT,PTEN_NULL, setdiff(as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL))) )
modelWT = lmerTest::lmer(data=dr, MODEL)

model = lmerTest::lmer(data=dr, MODEL)
model.rnef = lme4::ranef(model)


pvals.WT = coefficients(summary(modelWT))[,5]
pvals.GFP = coefficients(summary(modelGFP))[,5]
pvals = merge(data.frame(Genotype=names(pvals.WT), pvalue.from.WT=pvals.WT), data.frame(Genotype=names(pvals.GFP), pvalue.from.NULL=pvals.GFP), all.x=T, all.y=T)
###

## Add coefficients and CI
pvals <- merge(pvals, data.frame(Genotype=names(fixef(model)), Coefficients=fixef(model)), all.x=T, all.y=F)
model.confint =  data.frame( model.ci=confint(model) )
colnames(model.confint) <- paste0("axon",c("95CI.lo", "95CI.high"))
model.confint$Genotype = rownames(model.confint)

pvals <- merge(pvals, model.confint[grepl(rownames(model.confint), pattern = "Genotype"),], by="Genotype", all.x=T, all.y=T)
pvals$Genotype <- gsub(pvals$Genotype, pattern = "Genotype", replacement = "")

adjusted.means = adjust.data.for.batch.lme4(dr, model, 
                                            model.rnef = model.rnef, rnef.names = c("plate"),
                                            filter.agg = "axon", value.name="axon", pattern = "Genotype")
adjusted.means <- merge(adjusted.means, pvals, by="Genotype", all.x=T, all.y=F)

adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "axon", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)
adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "axon.adj", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)

writeObject(adjusted.means, filename = "axon.MEM.adjusted.means.df", row.names = F)
