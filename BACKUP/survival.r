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
	
  survival.eset = cbind(survival.data[colnames(exprs.data), c(time.col, event.col)], t(exprs.data))
    
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

