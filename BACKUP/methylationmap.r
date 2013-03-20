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
