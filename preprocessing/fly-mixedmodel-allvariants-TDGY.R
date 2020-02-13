# description: "Fly data mixed-effect model/preprocessing from the raw data with all variants in one model."
# author: "Manuel Belmadani"
# date: "21/02/2019"
source('utils.R')
MODEL = as.formula("time ~ Genotype +  (1 | vial)")

fly.indentity.raw <- read.csv(sep="\t",
                              here::here("Data/PTEN-Fly.identity.tsv"), 
                              header=1)

fly.indentity.raw$eclosion <- factor(fly.indentity.raw$eclosion)

dr <- as.tibble(fly.indentity.raw)
names(dr)[2] = "Genotype"

PTEN_WT = "WT"
PTEN_NULL = "attp2"
PTEN_INTERNAL = "TDGY"
dr$vial <- paste0(dr$Genotype,"_",dr$eclosion,"_",dr$vial)
dr$vial <- factor(dr$vial)

dr$Genotype <- as.character(dr$Genotype)
dr$Genotype[dr$condition == "Control"] <- PTEN_INTERNAL
dr$Genotype <- as.factor(dr$Genotype)

dr$WT<- dr$Genotype == PTEN_WT
dr$PTEN_NULL <- dr$Genotype == PTEN_NULL
dr$Control <- dr$PTEN_NULL | dr$WT

dr$Genotype <- factor(dr$Genotype,  c(setdiff(as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL)), PTEN_NULL, PTEN_WT) )
dr$condition <- factor(dr$condition,  c("Experiment","Control") )

PTEN_WT.pvalue = PTEN_WT
PTEN_NULL.pvalue = "conditionControl" # Tell is to use the conditionExperiment p-value for the Null

dr$Genotype = droplevels(dr$Genotype)

#Compute p-values from WT and NULL comparisons
dr$Genotype <- factor(dr$Genotype,
                      c( PTEN_NULL, PTEN_WT, setdiff(as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL))) )
modelATTP2 = lmerTest::lmer(data=dr, MODEL)

dr$Genotype <- factor(dr$Genotype,
                      c( PTEN_WT,PTEN_NULL, setdiff(as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL))) )
modelWT = lmerTest::lmer(data=dr, MODEL)

# Compute "standard model"
model = lmerTest::lmer(data=dr, MODEL)
model.rnef = lme4::ranef(model)

pvals.WT = coefficients(summary(modelWT))[,5]
modelATTP2 = coefficients(summary(modelATTP2))[,5]
pvals = merge(data.frame(Genotype=names(pvals.WT), pvalue.from.WT=pvals.WT), data.frame(Genotype=names(modelATTP2), pvalue.from.NULL=modelATTP2), all.x=T, all.y=T)


## Add coefficients and CI
pvals <- merge(pvals, data.frame(Genotype=names(fixef(model)), Coefficients=fixef(model)), all.x=T, all.y=F)
# model.confint =  data.frame( model.ci=confint(model) )
# colnames(model.confint) <- paste0("time",c("95CI.lo", "95CI.high"))
# model.confint$Genotype = rownames(model.confint)
# pvals <- merge(pvals, model.confint[grepl(rownames(model.confint), pattern = "Genotype"),], by="Genotype", all.x=T, all.y=T)

pvals$Genotype <- gsub(pvals$Genotype, pattern = "Genotype", replacement = "")

adjusted.means = adjust.data.for.batch.lme4(dr, model, 
                                            model.rnef = model.rnef, rnef.names = c("vial"),
                                            filter.agg = "time", value.name="time", pattern = "Genotype")
adjusted.means <- merge(adjusted.means, pvals, by="Genotype", all.x=T, all.y=F)

adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "time", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)
adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "time.adj", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)

writeObject(adjusted.means, filename = "fly.MEM.TDGY.adjusted.means.df", row.names = F)

