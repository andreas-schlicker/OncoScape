# Fits a COX model for each gene in the gene list.
# exprs.data should be a matrix with genes in the rows and samples in the columns
# survival.data should be matrix with samples in the rows and the annotation in columns
# genes should be a list of genes; the IDs used have to map to row names in the exprs.data matrix
# time.col should be the index or name of the column containing the survival time
# event.col should be the index or name of the column containing the event
coxphSurvival = function(exprs.data, survival.data, genes, time.col="dfs_time", event.col="dfs_event") {  
  if (!require(survival)) {
    stop("Can't load required package \"survival\".")
  }
	
  survival.eset = cbind(survival.data[colnames(exprs.data), c(time.col, event.col)], t(exprs.data[intersect(rownames(exprs.data), genes), ]))
    
  # Perform Cox regression for each gene and get the p-values
  cox.pvalues = matrix(NA, nrow=length(genes), ncol=2)
  colnames(cox.pvalues) = c("p.value", "coefficient")
  rownames(cox.pvalues) = genes
  for (i in genes) {
    # Try to calculate Cox regression; return NA if any error occured
    assign("coxMod", tryCatch(summary(coxph(Surv(get(time.col), get(event.col)) ~ get(i), data=survival.eset)), error = function(e) NA))
    if (class(coxMod) == "summary.coxph") {
      cox.pvalues[i,] = c(coxMod$waldtest["pvalue"], coxMod$coefficients[, 1])
    }
  }
  
  # Perform Benjamini-Hochberg multiple testing correction
  cox.pvalues = cbind(cox.pvalues, BH=p.adjust(cox.pvalues[, 1], method="BH"))
  
  return(cox.pvalues)
}

# Summarizes the results of a survival analysis. A gene needs to be significant
# in only one of the survival comparisons to receive a score of 1.;
# surv.res a matrix containing the coefficients and (preferably multiple testing
# corrected) p-values resulting from a Cox regression for each gene; genes are
# expected to be in the rows and coefficients and p-values in columns.
# genes a character vector containing the ids of the genes to analyze; ids have
# to map to row names in the "surv.res" matrix
# sig.cutoff the significance cut-off to be applied; defaults to 0.05
# direction either "down" or "up" defines whether negative or positive, 
# respectively, coefficients are expected
summarizeSurvival = function(surv.res, genes, sig.cutoff=0.05, direction=c("down", "up")) {
  direction = match.arg(direction)
 
  # Get the correct comparison function
  # If we want to find genes positively correlated with hazard, get the greaterThan function
  # If we want to find genes negatively correlated with hazard, get the smallerThan function
  compare = switch(direction, down=smallerThan, up=greaterThan)	
	
  gene.scores = rep(NA, length(genes))
  names(gene.scores) = genes
  
  steps = seq(1, ncol(surv.res), by=2)
  
  for (gene in genes) {
    score = 0
    for (i in steps) {
      if (!is.na(surv.res[gene, i]) && compare(surv.res[gene, i], 0) && !is.na(surv.res[gene, i+1]) && surv.res[gene, i+1] < sig.cutoff) {
        score = 1
      }
    }
    gene.scores[gene] = score
  }
  
  return(gene.scores)
}

