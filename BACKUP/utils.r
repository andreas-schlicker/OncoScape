# Perform Wilcoxon tests on the given matrix.
# If the matchedSamples argument is defined, a paired test is performed. In this
# case, the first length(matchedSamples) number of columns need to contain group
# 1 and then all matched samples in group 2 should be given. Samples in groups 1
# and 2 need to be in matched order.
# If a non-paired test is to be performed, groups should be a named vector of
# 1 and 2 indicating which samples belong to group 1 and 2.
doWilcox = function(inpMat, matchedSamples=NULL, groups=NULL) {
  # Paired Wilcoxon test?
  paired = !is.null(matchedSamples)
 
  # p-values from tests
  wilcox.p = rep(-1, times=nrow(inpMat))
  if (!is.null(rownames(inpMat)))
    names(wilcox.p) = rownames(inpMat)
 
  # Run paired Wilcoxon tests
  if (paired) {
    lms = length(matchedSamples)
    ncolMat = ncol(inpMat)
    for (i in 1:nrow(inpMat)) {
      wilcox.p[i] = tryCatch(wilcox.test(x=inpMat[i, 1:lms], y=inpMat[i, (lms+1):ncolMat], paired=TRUE)$p.value, error = function(e) NA)
    }
  } else {
    group1 = names(groups[groups == 1])
    group2 = names(groups[groups == 2])
    for (i in 1:nrow(inpMat)) {
      wilcox.p[i] = tryCatch(wilcox.test(x=inpMat[i, group1], y=inpMat[i, group2], paired=FALSE)$p.value, error = function(e) NA)
    }
  }
  return(wilcox.p)
}

# Calculate prioritization score by scoring each single dataset
# x should be a matrix with 10 columns
categoryScore = function(x) {
  # Final score
  res = 0
  # At least one Illumina 27 is higher methylated in tumors than normals
  # (significant in paired test) and negatively correlated with expression
  if (!is.na(x[1]) && x[1] > 0) {
    res = res + 1
  }
  # At least one Illumina 27 is higher methylated in tumors than normals
  # (significant in unpaired test) and negatively correlated with expression
  if (!is.na(x[2]) && x[2] > 0) {
    res = res + 1
  }
  # At least one Illumina 450 is higher methylated in tumors than normals
  # (significant in paired test) and negatively correlated with expression
  if (!is.na(x[3]) && x[3] > 0) {
    res = res + 1
  }
  # At least one Illumina 450 is higher methylated in tumors than normals
  # (significant in unpaired test) and negatively correlated with expression
  if (!is.na(x[4]) && x[4] > 0) {
    res = res + 1
  }
  # Copy number is significantly lower in tumors than in paired normals
  if (!is.na(x[5]) && x[5] < 0) {
    res = res + 1
  }
  # Gene is significantly associated with poorer disease free survival and
  # higher expression -> shorter survival
  if (!is.na(x[7]) && x[7] < 0.05 && !is.na(x[8]) && x[8] < 0) {
    res = res + 1
  }
  # Gene is significantly associated with poorer disease specific survival and
  # higher expression -> shorter survival
  if (!is.na(x[9]) && x[9] < 0.05 && !is.na(x[10]) && x[10] < 0) {
    res = res + 1
  }
  return(res)
}

# Calculate prioritization score by scoring each data type
# x should be a matrix with 10 columns
# sig.cutoff is the significance cut-off used
# direction either "down" or "up" for scoring evidence for the respective direction
# of regulation in tumors 
dataTypeScore = function(x, sig.cutoff=0.05, direction=c("down", "up")) {
  direction = match.arg(direction)
 
  # Get the correct comparison function
  # If we want to find genes down-regulated in tumor, get the greaterThan function
  # If we want to find genes up-regulated in tumor, get the smallerThan function
  compare = switch(direction, down=smallerThan, up=greaterThan)
  
  # Final score
  res = 0
  # At least one Illumina 27 probe is higher methylated in tumors than normals, 
  # significant in paired or unpaired test and negatively correlated with expression   
  # At least one Illumina 450 probe is higher methylated in tumors than normals, 
  # significant in paired or unpaired test and negatively correlated with expression
  if ((!is.na(x[1]) && x[1] > 0) || (!is.na(x[2]) && x[2] > 0) ||
      (!is.na(x[3]) && x[3] > 0) || (!is.na(x[4]) && x[4] > 0)) {
    res = res + 1
  }
  # Copy number is significantly lower in tumors than in paired normals
  if (!is.na(x[7]) && x[7] == 1) {
    res = res + 1
  }
  # Gene is significantly associated with poorer disease free or disease specific 
  # survival and lower expression -> shorter survival
  # Negative coefficients indicate lower expression is associated with shorter survival
  # Positive coefficients indicate higher expression is associated with shorter survival
  if ((!is.na(x[8]) && compare(x[8], 0) && !is.na(x[9]) && x[9] < sig.cutoff) || 
      (!is.na(x[10]) && compare(x[10], 0) && !is.na(x[11]) && x[11] < sig.cutoff)) {
    res = res + 1
  }
  return(res)
}

# Calculates the ratio of probes per gene that are significantly negatively associated with expression and all probes for that gene.
# meth.wilcox is a vector with p-values (preferably multiple-testing corrected) for methylation difference between tumor and normal
# meth.diff is a vector giving the value mean(methylation.tumor) - mean(methylation.normal) for each probe
# cors is a vector containing the correlation between methylation probes and expression of the corresponding gene
# meth.ann is a matrix mapping methylation probes to genes
# genes is a character vector of genes to consider
# wilcox.cutoff is the significance cut-off, probes with meth.wilcox < wilcox.cutoff will be considered significant
# diff.cutoff probes need to have a methylation difference smaller or greater than this cut-off to be considered significant
# regulation either "down" or "up" for finding genes that are regulated in the corresponding direction
# cor.min probes have to have a correlation larger or equal than this cut-off to be considered significant
# cor.max probes have to have a correlation smaller or equal than this cut-off to be considered significant
summarizeMethExprs = function(meth.wilcox, meth.diff, cors, meth.ann, genes, wilcox.cutoff=0.05, diff.cutoff=0.1, regulation=c("down", "up"), cor.min=-1.0, cor.max=-0.1) { 
  regulation = match.arg(regulation)
  
  # Get the correct comparison function
  # If we want to find genes with higher methylation in tumors get the greaterThan function
  # If we want to find genes with lower methylation in tumors, get the smallerThan function
  compare = switch(regulation, down=greaterThan, up=smallerThan)
  
  gene.scores = rep(0.0, length(genes))
  names(gene.scores) = genes
  
  for (gene in genes) {
    # Get all probes for that gene
    geneProbes = which(meth.ann[,1] == gene)
    if (length(geneProbes) > 1) {
      # If there is more than one probe, get the rownames 
      allProbes = rownames(meth.ann[geneProbes, ])
    } else {
      # Else the probe name has to be created again
      allProbes = c(paste(meth.ann[geneProbes, "chrom"], meth.ann[geneProbes, "pos"], sep="_"))
    }
    
    # Get all probes that show significantly different methylation levels between tumors and normals
    signifProbes = intersect(names(meth.wilcox)[meth.wilcox < wilcox.cutoff], allProbes)
    # Check which of these probes have a negative correlation with the expression data
    signifProbes = signifProbes[cors[signifProbes] <= cor.max & cors[signifProbes] >= cor.min]
    # Get the actual methylation differences
    diffs = meth.diff[signifProbes]
    # Get the ratio of (#significant probes with higher methylation) / (#all probes)
    gene.scores[gene] = length(diffs[compare(diffs, diff.cutoff)]) / length(allProbes) 
  }
  
  return(gene.scores)
}

# For each gene, summarizes evidence for significant copy number alterations and the correlation with gene expression. 
# cnv.wilcox is a vector with p-values (preferably multiple-testing corrected) for copy number difference between tumor and normal
# cnv.diff is a vector giving the value mean(cnv.tumor) - mean(cnv.normal) for each gene
# cors is a vector containing the correlation between copy number and expression for each gene
# genes is a character vector of genes to consider
# wilcox.cutoff is the significance cut-off, probes with cnv.wilcox > wilcox.cutoff will be considered significant
# diff.cutoff genes need to have a copy number difference smaller or greater than this cut-off to be considered significant
# regulation either "down" or "up" for finding genes that are lost or gained, respectively
# cor.min genes need to have a correlation larger or equal than this cut-off to be considered significant
# cor.max genes need to have a correlation smaller or equal than this cut-off to be considered significant
summarizeCnvExprs = function(cnv.wilcox, cnv.diff, cors, genes, wilcox.cutoff=0.05, diff.cutoff=0.1, regulation=c("down", "up"), cor.min=0.1, cor.max=1.0) { 
  regulation = match.arg(regulation)
  
  # Get the correct comparison function
  # If we want to find genes with higher copy number in tumors get the greaterThan function
  # If we want to find genes with lower copy number in tumors, get the smallerThan function
  compare = switch(regulation, down=smallerThan, up=greaterThan)
  
  gene.scores = rep(0.0, length(genes))
  names(gene.scores) = genes
  
  for (gene in genes) {
    # Get p-value for the copy number difference for that gene
    significant = cnv.wilcox[gene] < wilcox.cutoff
    # Check whether the copy number is correlated with expression
    correlated = (cors[gene] >= cor.min && cors[gene] <= cor.max)
    # Get the actual methylation differences
    difference = compare(cnv.diff[gene], diff.cutoff)
    # Does the gene pass all filters?
    gene.scores[gene] = as.integer(significant && correlated && difference) 
  }
  
  return(gene.scores)
}

greaterThan = function(x, y) { x > y }
smallerThan = function(x, y) { x < y }

