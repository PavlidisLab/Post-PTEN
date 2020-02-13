
DEFAULT.LABELS = list(
  Fly.activity = "Fly activity",
  Yeast.activity = "Yeast activity",
  Worm.activity = "Worm activity",
  PSD.95.Density = "PSD-95 density",
  Total.Dendrite.Length = "Total dendrite length",
  Axon.activity = "DRG axon outgrowth",
  
  Yeast.stability = "Yeast stability",
  Fly.stability = "Fly stability",
  Hek.stability = "HEK293 stability", 
  
  mighell_bin = "Mighell et al. (Activity)", 
  AssayMean = "Averaged activity", 
  matreyek_score = "Matreyek et al. (Stability)",
  CADD13Raw = "CADD Score (1.3)",
  SNAP2 = "SNAP2 Score",
  
  Hek.activity="HEK293 pAKT activity",
  Size.adj.scaled.means="Yeast activity",
  time.adj.scaled.means="Fly activity",
  Chemotaxis_Index.adj.scaled.means="Worm activity",
  PSD.95.Density.adj.scaled.means="PSD-95 Density",
  Total.Dendrite.Length.adj.scaled.means="Total dendrite length",
  axon.adj.scaled.means="DRG axon outgrowth"
  
)

SHORT.LABELS = list(
  Fly.activity = "Fly",
  Yeast.activity = "Yeast",
  Worm.activity = "Worm",
  PSD.95.Density = "PSD-95",
  Total.Dendrite.Length = "TDL",
  Axon.activity = "Axon",
  
  Yeast.stability = "Yeast (WB)",
  Fly.stability = "Fly (WB)",
  Hek.stability = "HEK293",
  
  Hek.activity="pAKT",
  Size.adj.scaled.means="Yeast",
  time.adj.scaled.means="Fly",
  Chemotaxis_Index.adj.scaled.means="Worm",
  PSD.95.Density.adj.scaled.means="PSD-95",
  Total.Dendrite.Length.adj.scaled.means="TDL",
  axon.adj.scaled.means="Axon"
)

VERBOSE.LABELS = list(
  Fly.activity = "Fly eclosion",
  Yeast.activity = "Yeast growth",
  Worm.activity = "Worm chemotaxis",
  PSD.95.Density = "PSD-95 synapse density",
  Total.Dendrite.Length = "Total dendrite length",
  Axon.activity = "DRG axon outgrowth",
  
  Yeast.stability = "Yeast western blot",
  Fly.stability = "Fly western blot",
  Hek.stability = "HEK293 protein abundance",
  
  Hek.activity="HEK293 Phospho-AKT",
  Size.adj.scaled.means="Yeast growth",
  time.adj.scaled.means="Fly eclosion",
  Chemotaxis_Index.adj.scaled.means="Worm chemotaxis",
  PSD.95.Density.adj.scaled.means="PSD-95 synapse density",
  Total.Dendrite.Length.adj.scaled.means="Total dendrite length",
  axon.adj.scaled.means="DRG axon outgrowth"
)


CURRENT.LABELS = DEFAULT.LABELS
 
as.label <- 
  # Simply translate using CURRENT.LABELS 
  
  # Examples:
  # as.label(c("Yeast.activity", "Fly.activity"))
  # as.label(colnames(combined))

  function(V, MAP=CURRENT.LABELS){
    bool.factor = F
    order.factor = NULL
    if (is.factor(V)){
      bool.factor = T
      order.factor = as.numeric(V)
      V = as.character(V)
    }
    
    V.new = as.character(sapply(V, FUN = function(x) {
      if (x %in% names(MAP)) {
        return(MAP[[x]]) 
      } else {
        return(x)
      }
    }))
    
    if (bool.factor){
      V.new = as.factor(V.new)
      V.new = fct_reorder(V.new, order.factor)
    }
    
    return(V.new)
}


