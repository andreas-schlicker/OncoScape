# Compares methylation data of tumor and normal samples.
# The function will run paired and un-paired Wilcoxon tests. Pairing of samples
# is determined by matching sample names in the tumors and normals matrices.
# tumors is a methylation data matrix with samples in the columns and probes in the rows
# tumors.probeann is a matrix with annotation of the probes;
# either probes are in the same order as in tumors or the rownames match
# normals is a methylation data matrix with samples in the columns and probes in the rows
# genes is a vector of genes to be tested; gene ids have to match the first column of the probe annotation matrix
# Returns a named list with results from paired and un-paired tests.
runMethComp = function(tumors, tumors.probeann, normals, genes) {
  # Samples from tumors that have a matched normal
  tumors.matchedsamples = intersect(colnames(tumors), colnames(normals))
  # Rename normal samples to reflect the state
  colnames(normals) = paste(colnames(normals), "normal", sep="_")

  tumors.groups = c(rep(1, times=ncol(tumors)), rep(2, times=ncol(normals)))
  names(tumors.groups) = c(colnames(tumors), colnames(normals))

  # Which probes map to any of the CIS genes?
  selected.probes = rownames(tumors.probeann[which(tumors.probeann[,1] %in% genes), ])
  selected.probes.selected = intersect(intersect(rownames(tumors), rownames(normals)), selected.probes)

  # Run paired and un-paired Wilcoxon tests
  pairedwilcox.p = doWilcox(inpMat=cbind(tumors[selected.probes.selected, tumors.matchedsamples], normals[selected.probes.selected, ]), matchedSamples=tumors.matchedsamples)
  wilcox.p = doWilcox(inpMat=cbind(tumors[selected.probes.selected, ], normals[selected.probes.selected, ]), groups=tumors.groups)

  return(list(unpaired=wilcox.p, paired=pairedwilcox.p))
}

# Impute missing values 
# This method will impute missing data using the impute.knn function as implemented 
# in the "impute" library. If this library is not installed or impute=FALSE, all
# probes with missing values will be removed.
# meth.data is the matrix will methylation probes in rows and samples in columns
# impute is a boolean indicating whether missing values should be imputed
# no.na is the threshold giving the number of missing values from which on a
# probe will be removed. By default, a probe is only removed if its value is
# missing for all samples.
cleanMethylation = function(meth.data, impute=TRUE, no.na=(ncol(meth.data)-1)) {
  # How many values are missing for each probe?
  no.nas = apply(meth.data, 1, function(x) { sum(is.na(x)) })
  if (!require(impute) || !impute) {
    print("Could not load library \"impute\". All probes with missing values will be removed")
    # Remove probes with missing values
    exclude = which(no.nas > 0)
    meth.data.imputed = list(data=meth.data[-exclude, ])
  } else {
    # Probes to be excluded
    exclude = which(no.nas > no.na)
    # Impute missing values
    meth.data.imputed = impute.knn(meth.data[-exclude, ])
  }
  
  return(meth.data.imputed)
}

# Calculates correlation between methylation probes and expression of the corresponding genes.
# Correlation is calculated using all samples contained in both data matrices
# meth.data should be a matrix with probes in rows and samples in columns
# meth.ann has to be an annotation matrix mapping probes to genes; gene ids have to match
# ids used in the expression data.
# exprs.data should be a matrix with probes in rows and samples in columns
# meth.probes should be a character vector giving the probes to test
corMethExprs = function(meth.data, meth.ann, exprs.data, meth.probes) {
  # All genes contained in the expression data
  exprs.genes = rownames(exprs.data)
  commonSamples = intersect(colnames(meth.data), colnames(exprs.data))
	
  cors = c()
  n = c()
  for (x in meth.probes) {
    if (meth.ann[x, 1] %in% exprs.genes) {
     cors = c(cors, cor(exprs.data[meth.ann[x, 1], commonSamples], meth.data[x, commonSamples], method="spearman"))
     n = c(n, x)
    }
  }
  names(cors) = n
  
  return(cors)
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

# Summarizes the results of analyzing methylation data.
# methylation matrix containing the results of at least one methylation set. 
# Typically this will be the output of the "summarizeMethExprs" function.
# A gene receives a score of one if at least on entry for that gene exceeds the
# given threshold.
# genes a character vector containing the ids of the genes to analyze; ids have
# to map to row names in the "methylation" matrix
# threshold floating point number giving the threshold that has to be met;
# defaults to 0
# score either "atleast" (default) or "absolute". If set to "atleast", a gene
# receives a score of 1 if at least one column exceeds the threshold. If set to
# "absolute", the score of the gene will be equal to the number of columns 
# exceeding the threshold.
summarizeMethylation = function(methylation, genes, threshold=0, score=c("atleast", "absolute")) {
  score = match.arg(score)
  score = switch(score, atleast=TRUE, absolute=FALSE)
	
  # gene.scores = rep(0, length(genes))
  # names(gene.scores) = genes
  # 
  # nocol = ncol(methylation)
  # 
  # for (gene in genes) {
  #   score = 0
  #   for (i in nocol) {
  #     if (!is.na(methylation[gene, i]) && methylation[gene, i] > 0) {
  #       score = 1
  #     }
  #   }
  #   genes.score[gene] = score
  # }
  
  # How many comparisons per gene exceed the threshold? 
  gene.scores = apply(methylation, 1, function(x) { length(which(x > threshold)) })	
  # Is there at least one comparison exceeding the threshold?
  if (score)
    gene.scores = sapply(gene.scores, function(x) { as.integer(x > 0) }) 
  
  return(gene.scores)
}

