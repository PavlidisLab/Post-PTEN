# description: "Yeast data linear (model/preprocessing from the raw data."
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
source('utils.R')
MODEL = as.formula("Size ~ Genotype + (1|paired)")
# MODEL = as.formula("Size ~ Genotype + (1|paired)")
# MODEL = as.formula("Size ~ Genotype  + (1|paired) ")


PTEN_WT = "WT" 
PTEN_NULL = "pEGH"

yeast.raw<-read.csv(sep="\t",
                    here::here("Data/PTEN-MEM-Yeast.tsv"),
                    header=1)

yeast.raw <- yeast.raw[yeast.raw$Gene.Name==YEAST.SENTINEL, 
                       c("Variant", "VariantStatus", "Screen",  "Size", "paired")]

yeast.raw$Variant = gsub(pattern = "WT.*", replacement = PTEN_WT, x = yeast.raw$Variant)
yeast.raw$Variant = gsub(pattern = "pEGH.*", replacement = PTEN_NULL, x = yeast.raw$Variant)
yeast.raw$Variant[yeast.raw$VariantStatus == "ControlSpot"] <- PTEN_WT

yeast.raw$PTENVariant <- factor(yeast.raw$Variant)
yeast.raw$Screen <- factor(yeast.raw$Screen)
yeast.raw$paired <- factor(yeast.raw$paired)

dr <- as.tibble(yeast.raw)
dr
names(dr)[1] = "Genotype"

dr$WT<- dr$Genotype == "WT"
dr$PTENNull <- grepl(pattern = PTEN_NULL, dr$Genotype)
dr$Control <- dr$PTENNull | dr$WT
dr$Genotype <- factor(dr$Genotype)
yeast.adjusted.means.df = NULL
dr$Genotype = droplevels(dr$Genotype)
dr$RawSize = dr$Size

########################################################
for ( plate in unique(dr$Screen) ){
  # dr[dr$Screen == plate, "Size"] <- scale(dr[dr$Screen == plate, "Size"])
  # dr[dr$Screen == plate, "Size"] <- dr[dr$Screen == plate, "Size"] / median(unlist( dr[dr$Screen == plate & dr$Genotype == "pEGH", "Size"]) )
  # dr[dr$Screen == plate, "Size"] <- dr[dr$Screen == plate, "Size"] - mean(unlist( dr[dr$Screen == plate & dr$Genotype == "WT", "Size"]) )
    # dr[dr$Screen == plate, "Size"] <- dr[dr$Screen == plate, "Size"] - mean(unlist( dr[dr$Screen == plate & dr$VariantStatus == "ControlSpot", "Size"]) )
  
  # mean.wt = median(dr$Size[dr$Screen == plate & dr$Genoype == "WT"])
  # dr$Size[plate$Screen == plate] <- dr$Size[yeast.multiplate$Screen == plate]  - mean.wt

  mean.ev = median(dr$Size[dr$Screen == plate & dr$Genotype == "pEGH"])
  dr$Size[dr$Screen == plate] <- dr$Size[dr$Screen == plate]  / mean.ev
}
#######################################################


# Apply mixed effect model; treat plate as a mixed effect
# Compute p-values from WT and NULL comparisons
dr$Genotype <- factor(dr$Genotype,
                      c(PTEN_NULL, PTEN_WT, setdiff(
                        as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL)
                      )))
#######################
# dr = backup.dr
# # C124S.4A.plate = unique( dr$Screen[dr$Genotype == "G44D"] )
# C124S.4A.plate = unique( dr$Screen[dr$Genotype == "C124S-4A"] )
# # dr = subset(dr, dr$Screen == C124S.4A.plate)
# # dr = subset(dr, dr$Genotype %in%  c("pEGH","WT","G44D"))
# dr = subset(dr, dr$Genotype %in%  c("pEGH","WT","C124S-4A"))
# dr = dr[dr$paired %in% dr$paired[duplicated(dr$paired)],]
# dr = dr[ !(dr$WT & dr$VariantStatus != "ControlSpot"),] # Remove non controlspot WTs
# summary(dr)
# summary(dr$Genotype)

######################

modelNULL = lmerTest::lmer(data=dr, MODEL)
# modelNULL = robustlmm::rlmerRcpp(data=dr, MODEL)

dr$Genotype <- factor(dr$Genotype,
                      c(PTEN_WT, PTEN_NULL, setdiff(
                        as.character(dr$Genotype), c(PTEN_WT, PTEN_NULL)
                      )))
modelWT = lmerTest::lmer(data=dr, MODEL)
# modelWT = robustlmm::rlmerRcpp(data=dr, MODEL)

model = lmerTest::lmer(data=dr, MODEL)
# model = robustlmm::rlmerRcpp(data=dr, MODEL)
model.rnef = lme4::ranef(model)

pvals.WT = coefficients(summary(modelWT))[,5]
pvals.NULL = coefficients(summary(modelNULL))[,5]
pvals = merge(data.frame(Genotype=names(pvals.WT), pvalue.from.WT=pvals.WT),
              data.frame(Genotype=names(pvals.NULL), pvalue.from.NULL=pvals.NULL),
              all.x=T, all.y=T)
###

pvals <- merge(pvals, data.frame(Genotype=names(fixef(model)), Coefficients=fixef(model)), all.x=T, all.y=F) # Add coefficients and CI

model.confint =  data.frame( model.ci=confint(model) )
colnames(model.confint) <- paste0("Size",c("95CI.lo", "95CI.high"))
model.confint$Genotype = rownames(model.confint)
pvals <- merge(pvals, model.confint[grepl(rownames(model.confint), pattern = "Genotype"),], by="Genotype", all.x=T, all.y=T)

pvals$Genotype <- gsub(pvals$Genotype, pattern = "Genotype", replacement = "")

adjusted.means = adjust.data.for.batch.lme4(dr, model, 
                                            model.rnef = model.rnef, 
                                            rnef.names = c("paired"),
                                            value.name="Size", 
                                            pattern = "Genotype")

adjusted.means <- merge(adjusted.means, pvals, by="Genotype", all.x=T, all.y=F)

adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "Size", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)
adjusted.means = add.scaled.01.means(df = adjusted.means, measurement = "Size.adj", ctrl.null = PTEN_NULL, ctrl.wt = PTEN_WT)
#################################################################################################################################
writeObject(adjusted.means, filename = paste0("yeast.MEM.",YEAST.SENTINEL,".adjusted.means.df"), row.names = F)


