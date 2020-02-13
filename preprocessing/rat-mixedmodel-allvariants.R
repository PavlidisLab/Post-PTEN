# description: "Preprocessing the Raw Rat data using a linear mixed-effects model from the raw data with all variant in one model."
# author: "Manuel Belmadani"
# date: "23/04/2018"
# updated: "23/04/2018"
library("lme4")

source("utils.R")

rat.data <- read.csv(sep = "\t", 
                     file = "Data/PTEN-Rat.tsv")

rat.data[rat.data=="TBC"] <-NA
rat.data[rat.data=="NA"] <-NA
rat.data[rat.data==""] <-NA
rat.data[,4:7] <- apply(MARGIN = 2, X=rat.data[,4:7], FUN=as.numeric)
rat.data = cbind(rat.data, PSD.95.Density.log=log10(rat.data$PSD.95.Density), Total.Dendrite.Length.log10=log10(rat.data$Total.Dendrite.Length))

rat.raw.df <- rat.data
summary(rat.raw.df)

rat.pvals.df <- NULL
rat.coef.df <- NULL
rat.stderr.df <- NULL
rat.adjusted.means.df = NULL

for ( i in c(4,5,6,7,8,9)){
  dr <- as.tibble(rat.raw.df)
  measure.name <- colnames(rat.raw.df)[i]

  PTEN_WT = "Hum_WT"
  PTEN_NULL = "GFP_Control"
  
  dr$PTEN_NULL<- dr$Condition == PTEN_NULL
  dr$WT<- dr$Condition == PTEN_WT
  
  dr$Control<- dr$WT | dr$PTEN_NULL
  dr$Condition <- factor(dr$Condition,  c(setdiff(as.character(dr$Condition), c(PTEN_WT, PTEN_NULL)), PTEN_NULL, PTEN_WT) )
  
  names(dr)[i] <- "Measure"
  
  x = dr
  x <- as.data.frame(x[!is.na(x$Measure),]) # Remove NAs
  x$Condition = droplevels(x$Condition)
  model = NULL
    
  if (length(unique(x$Culture)) > 1) {
    C = "  + ( 1 | Culture ) "
    
  } else {
    C = ""
  }
  
  if (length(unique(x$Measure_20X_Masking)) > 1) {
    M20X = " + ( 1 | Measure_20X_Masking ) "
  } else {
    M20X = ""
  }
  
  if (length(unique(x$Measure_60X_Masking)) > 1) {
    M60X = " + ( 1 | Measure_60X_Masking ) "
  } else {
    M60X = ""
  }
  
  if (length(unique(x$Measure_Confocal_Imaging)) > 1) {
    MCI = " + ( 1 | Measure_Confocal_Imaging ) "
  } else {
    MCI = ""
  }
  
  
  if (length(unique(x$Measure_PSD.95_Density)) > 1) {
    MPSD95 = " + ( 1 | Measure_PSD.95_Density ) "
  } else {
    MPSD95 = ""
  }
  
  if (length(unique(x$Measure_Gephyrin_Density)) > 1) {
    MG = " + ( 1 | Measure_Gephyrin_Density ) "
  } else {
    MG = ""
  }
  
  if (length(unique(x$Measure_Total_Dendrite_Length)) > 1) {
    MTDL = " + ( 1 | Measure_Total_Dendrite_Length ) "
  } else {
    MTDL = ""
  }
    
  if (i == 4 || i == 8 ) {
    variantFormula = paste0("Measure ~ Condition", C, M60X, MCI, MPSD95)
    variantFormula.null = paste0("Measure ~ 1", C, M60X, MCI, MPSD95)
  } else if (i == 5 || i == 9) {
    variantFormula = paste0("Measure ~ Condition", C, M60X, MCI, MG)
    variantFormula.null = paste0("Measure ~ 1", C, M60X, MCI, MG)
  } else if (i == 6) {
    variantFormula = paste0("Measure ~ Condition", C, M20X, MCI, MTDL)
    variantFormula.null = paste0("Measure ~ 1", C, M20X, MCI, MTDL)
  } else if (i == 7) {
    variantFormula = paste0("Measure ~ Condition", C, M20X, MCI)
    variantFormula.null = paste0("Measure ~ 1", C, M20X, MCI)
  }
    
  SIMPLE = T # Comparing the models with an anova revealed that the model fit is not better/possibly worse when including the experimenter.
  if (SIMPLE) {
    # Use simple version of fixed effect model.
    CULTURE = "  + (1|Culture) "
    MODEL = paste0("Measure ~ Condition", CULTURE)
    
    x$Condition <- factor(x$Condition,
                          c( PTEN_NULL, PTEN_WT, setdiff(as.character(x$Condition), c(PTEN_WT, PTEN_NULL))) )
    modelGFP = lmerTest::lmer(data=x, MODEL)
    x$Condition <- factor(x$Condition,
                          c( PTEN_WT,PTEN_NULL, setdiff(as.character(x$Condition), c(PTEN_WT, PTEN_NULL))) )
    modelWT = lmerTest::lmer(data=x, MODEL)
    
    model = lmerTest::lmer(data=x, MODEL)
    model.rnef = lme4::ranef(model)
    
    
    pvals.WT = coefficients(summary(modelWT))[,5]
    pvals.GFP = coefficients(summary(modelGFP))[,5]
    pvals = merge(data.frame(Condition=names(pvals.WT), pvalue.from.WT=pvals.WT), data.frame(Condition=names(pvals.GFP), pvalue.from.NULL=pvals.GFP), all.x=T, all.y=T)
    
    ## Add coefficients and CI
    pvals <- merge(pvals, data.frame(Condition=names(fixef(model)), Coefficients=fixef(model)), all.x=T, all.y=F)
    model.confint =  data.frame( model.ci=confint(model) )
    colnames(model.confint) <- paste0("Measure", c("95CI.lo", "95CI.high"))
    model.confint$Condition = rownames(model.confint)
    
    pvals <- merge(pvals, model.confint[grepl(rownames(model.confint), pattern = "Condition"),], by="Condition", all.x=T, all.y=T)
    pvals$Condition <- gsub(pvals$Condition, pattern = "Condition", replacement = "")
    
    adjusted.means = adjust.data.for.batch.lme4(
      x,
      model,
      value.name = "Measure",
      filter.agg = "Measure",
      pattern = "Condition",
      variant.name = "Condition",
      agg.by.list = T
    )
    adjusted.means = adjusted.means[, apply(adjusted.means, MARGIN = 2, function(x) sum(is.na(x)) != length(x))] # Drop NA cols
    adjusted.means = adjusted.means[, grepl(names(adjusted.means), pattern = "Condition|Measure.adj|^Measure$|^sd.Measure$")] # Drop other redundant cols
    
    adjusted.means = add.scaled.01.means(
      df = adjusted.means,
      measurement = "Measure.adj",
      ctrl.null = PTEN_NULL,
      ctrl.wt = PTEN_WT,
      select.by = "Condition"
    )
    
    adjusted.means <- merge(adjusted.means, pvals, by.x="Condition", by.y="Condition", all.x=T, all.y=F)
    names(adjusted.means)[names(adjusted.means) == "Condition"] = "Genotype" #TODO: Should be done earlier. No reason to name it condition.
    
    adjusted.means$Phenotype = measure.name
    if (is.null(rat.adjusted.means.df)) {
      rat.adjusted.means.df <<- adjusted.means
    } else {
      rat.adjusted.means.df <<-
        rbind(rat.adjusted.means.df, adjusted.means)
    }
    
  } else {
    # Use mixed effect model with experimenter information
    stop("Mixed effect model with experimenter information is disabled.")
    # if (as.logical(grepl("|", x = variantFormula, fixed = T))) {
    #   modelBy <-
    #     function(formula = NULL, data = NULL)
    #       lmer(data = data, formula, REML = F)
    # } else {
    #   modelBy <- lm
    # }
    # 
    # model <- modelBy(data = x, formula = as.formula(variantFormula))
    # model.null <-
    #   modelBy(data = x, formula = as.formula(variantFormula.null))
    # 
    # anova.model <-
    #   anova(model, model.null)  # Use x[8][2,1] for p-value
    # anova.model
    # 
    # model.ci = confint(model) # TODO: This fails and is pretty expensive; temporarily fallback to SIMPLE model.
    # model.effect = cbind(summary(model)$coefficient, pvalue = anova.model[8]$`Pr(>Chisq)`)
    # cbind(model.effect, model.ci[row.names(model.effect), ])
  }
}

writeObject(rat.adjusted.means.df, filename = "rat.MEM.adjusted.means.df", row.names = F)