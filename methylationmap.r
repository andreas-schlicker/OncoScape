# Compares methylation data of tumor and normal samples.
# The function will run paired and un-paired Wilcoxon tests. Pairing of samples
# is determined by matching sample names in the tumors and normals matrices.
# tumors is a methylation data matrix with samples in the columns and probes in the rows
# tumors.probeann is a matrix with annotation of the probes;
# either probes are in the same order as in tumors or the rownames match
# normals is a methylation data matrix with samples in the columns and probes in the rows
# genes is a vector of genes to be tested; gene ids have to match the first column of the probe annotation matrix
# Returns a named list with results from paired and un-paired tests.

runMethComp = function(tumors, normals, probes) {
  # Samples from tumors that have a matched normal
  tumors.matchedsamples = intersect(colnames(tumors), colnames(normals))
  # Retain only the normal samples for which there is a tumor
  normals = normals[, tumors.matchedsamples]
  # Rename normal samples to reflect the state
  colnames(normals) = paste(colnames(normals), "normal", sep="_")

  tumors.groups = c(rep(1, times=ncol(tumors)), rep(2, times=ncol(normals)))
  names(tumors.groups) = c(colnames(tumors), colnames(normals))

  # Which probes map to any of the genes?
  selected.probes = intersect(intersect(rownames(tumors), rownames(normals)), probes)

  inpMat = cbind(tumors[selected.probes, tumors.matchedsamples], normals[selected.probes, ])
  if (length(selected.probes) == 1) {
    inpMat = matrix(c(tumors[selected.probes, tumors.matchedsamples], normals[selected.probes, ]), nrow=1)
    rownames(inpMat) = selected.probes
    colnames(inpMat) = c(tumors.matchedsamples, colnames(normals))
  }
  # Run paired Wilcoxon tests
  pairedwilcox.p = doWilcox(inpMat=inpMat, matchedSamples=tumors.matchedsamples)
  
  inpMat = cbind(tumors[selected.probes, ], normals[selected.probes, ])
  if (length(selected.probes) == 1) {
    inpMat = matrix(c(tumors[selected.probes, ], normals[selected.probes, ]), nrow=1)
    rownames(inpMat) = selected.probes
    colnames(inpMat) = c(colnames(tumors), colnames(normals))
  }
  # Run un-paired Wilcoxon tests
  wilcox.p = doWilcox(inpMat=inpMat, groups=tumors.groups)

  return(list(unpaired=wilcox.p, paired=pairedwilcox.p))
}

# Impute missing values 
# This method will impute missing data using the impute.knn function as implemented 
# in the "impute" library. If this library is not installed or impute=FALSE, all
# probes with missing values will be removed.
# meth.data is the matrix with methylation probes in rows and samples in columns
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
	
  cors = numeric(nrow(meth.probes))
  n = character(nrow(meth.probes))
  for (i in 1:length(meth.probes)) {
    if (meth.ann[meth.probes[i], 1] %in% exprs.genes) {
      cors[i] = cor(exprs.data[meth.ann[meth.probes[i], 1], commonSamples], 
			  		meth.data[meth.probes[i], commonSamples], 
					method="spearman")
      n[i] = meth.probes[i]
    }
  }
  cors = cors[1:i]
  names(cors) = n[1:i]
  
  cors
}

# Calculates correlation between methylation probes and expression of the corresponding genes.
# This method deals with non-unique mappings of probes to genes by testing all
# of them.
# Correlation is calculated using all samples contained in both data matrices
# meth.data should be a matrix with probes in rows and samples in columns
# meth.ann has to be an annotation matrix mapping probes to genes;
# probe ids in 1st column, gene ids in 2nd;
# gene ids have to match ids used in the expression data.
# exprs.data should be a matrix with probes in rows and samples in columns
# meth.probes should be a character vector giving the probes to test
corMethExprsMult = function(meth.data, meth.ann, exprs.data, meth.probes, filt=c("all", "low", "high")) {
	filt = match.arg(filt)
	
	# All genes contained in the expression data
	exprs.genes = rownames(exprs.data)
	commonSamples = intersect(colnames(meth.data), colnames(exprs.data))
	samps = commonSamples
	
#	cors = c()
#	n = c()
#	for (x in meth.probes) {
#		if (filt == "low") {
#			samps = colnames(meth.data)[which(meth.data[x, commonSamples] < 0.3)]
#		} else if (filt == "high") {
#			samps = colnames(meth.data)[which(meth.data[x, commonSamples] > 0.7)]
#		}
#		genes = unlist(strsplit(meth.ann[which(meth.ann[, 1] == x), 2], ";"))
#		regions = unlist(strsplit(meth.ann[which(meth.ann[, 1] == x), 3], ";"))
#		genes_regions = unique(paste(genes, regions, sep="_"))
#		
#		for (gene_region in genes_regions) {
#			gs = strsplit(gene_region, "_")
#			if (gs[[1]][1] %in% exprs.genes) {
#				cors = c(cors, cor(exprs.data[gs[[1]][1], samps], meth.data[x, samps], method="spearman"))
#				n = c(n, paste(x, gene_region, sep="_"))
#			}
#		}
#	}
#	names(cors) = n
	
	# Correlations
	cors = list()
	# Names of correlations
	n = list()
	# Length of temporary vectors
	LENGTH = 1000
	# Temporary vectors
	temp1 = numeric(LENGTH)
	temp2 = character(LENGTH)
	# Counter for current element in temporary vector
	i = 1
	# Counter of temporary vectors
	j = 1
	
	for (x in meth.probes) {
		# No space left in temporary vectors
		# Save the vectors in the lists
		if (i > LENGTH) {
			cors[[j]] = temp1
			temp1 = numeric(LENGTH)
			n[[j]] = temp2
			temp2 = character(LENGTH)
			j = j + 1
			i = 1
		}
		
		# Potentially filter samples
		if (filt == "low") {
			samps = colnames(meth.data)[which(meth.data[x, commonSamples] < 0.3)]
		} else if (filt == "high") {
			samps = colnames(meth.data)[which(meth.data[x, commonSamples] > 0.7)]
		}
		# Get all genes and gene regions and combine them
		genes = unlist(strsplit(meth.ann[which(meth.ann[, 1] == x), 2], ";"))
		regions = unlist(strsplit(meth.ann[which(meth.ann[, 1] == x), 3], ";"))
		genes_regions = unique(paste(genes, regions, sep="_"))
		
		# Cycle through all gene-region combinations
		for (gene_region in genes_regions) {
			gs = strsplit(gene_region, "_")
			# If there is expression data, save the data in the temporary vectors
			if (gs[[1]][1] %in% exprs.genes) {
				temp1[i] = cor(exprs.data[gs[[1]][1], samps], meth.data[x, samps], method="spearman")
				temp2[i] = paste(x, gene_region, sep="_")
				i = i + 1
			}
		}
	}
	# Save the last values
	cors[[j]] = temp1[1:i]
	n[[j]] = temp2[1:i]
	# And build the return value
	cors = unlist(cors)
	names(cors) = unlist(n)
	
	cors
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

doMethylationAnalysis = function(tumors, 
								normals, 
								genes, 
								exprs, 
								probe.annotation, 
								samples=NULL, 
								wilcox.cutoff=0.05, 
								diff.cutoff=0.1, 
								regulation=c("down", "up"), 
								cor.min=0.1, 
								cor.max=1.0) {
							
  regulation = match.arg(regulation)
  
  # Get the correct comparison function
  # If we want to find genes down-regulated in tumors get the greaterThan function
  # If we want to find genes up-regulated in tumors, get the smallerThan function
  compare = switch(regulation, down=greaterThan, up=smallerThan)
  
  if (is.null(samples)) {
    tumor.samples = colnames(tumors)
    normal.samples = colnames(normals)
  } else {
    tumor.samples = intersect(samples, colnames(tumors))
    normal.samples = intersect(samples, colnames(normals))
  }
  
  selected.probes = intersect(intersect(rownames(tumors), rownames(normals)), rownames(probe.annotation[which(probe.annotation[,1] %in% genes), ]))
  
  mean.tumor = apply(tumors[selected.probes, tumor.samples], 1, mean)
  mean.normal = apply(normals[selected.probes, normal.samples], 1, mean)
  mean.diff = mean.tumor - mean.normal
  selected.probes = names(mean.diff[compare(mean.diff, diff.cutoff)])
  selected.probes = selected.probes[!is.na(selected.probes)]
  
  cors = c()
  if (length(selected.probes) > 0) {
    temp = matrix(tumors[selected.probes, tumor.samples], nrow=length(selected.probes), byrow=TRUE)
    rownames(temp) = selected.probes
    colnames(temp) = tumor.samples
    cors = corMethExprs(temp, probe.annotation, exprs, selected.probes)
    selected.probes = names(cors[which(cors <= cor.max & cors >= cor.min)])
  }
  
  wilcox = list(unpaired=c(), paired=c())
  wilcox.bh = list(unpaired=c(), paired=c())
  selected.probes.unpaired = c()
  selected.probes.paired = c()
  if (length(selected.probes) > 0) {
    tt = tumors[selected.probes, tumor.samples]
    tn = normals[selected.probes, normal.samples]
    if (length(selected.probes) == 1) {
      tt = matrix(tumors[selected.probes, tumor.samples], nrow=1, byrow=TRUE)
      colnames(tt) = tumor.samples
      rownames(tt) = selected.probes
      
      tn = matrix(normals[selected.probes, normal.samples], nrow=1, byrow=TRUE)
      colnames(tn) = normal.samples
      rownames(tn) = selected.probes
    }
  
    wilcox = runMethComp(tt, tn, selected.probes)
    wilcox.bh = list(unpaired=p.adjust(wilcox$unpaired, method="BH"), paired=p.adjust(wilcox$paired, method="BH"))
    selected.probes.unpaired = names(wilcox.bh$unpaired[which(wilcox.bh$unpaired < wilcox.cutoff)])
    selected.probes.paired = names(wilcox.bh$paired[which(wilcox.bh$paired < wilcox.cutoff)])
  }
  
  gene.scores = list(unpaired=rep(0.0, length(genes)), paired=rep(0.0, length(genes))) 
  names(gene.scores$unpaired) = genes
  names(gene.scores$paired) = genes
  
  for (gene in genes) {
    # Get all probes for that gene
    geneProbes = which(probe.annotation[,1] == gene)
    if (length(geneProbes) > 1) {
      # If there is more than one probe, get the rownames 
      allProbes = rownames(probe.annotation[geneProbes, ])
    } else {
      # Else the probe name has to be created again
      allProbes = c(paste(probe.annotation[geneProbes, "chrom"], probe.annotation[geneProbes, "pos"], sep="_"))
    }
    
    temp = length(geneProbes)
    if (temp > 0) {
      gene.scores$unpaired[gene] = length(intersect(allProbes, selected.probes.unpaired)) / temp
      gene.scores$paired[gene] = length(intersect(allProbes, selected.probes.paired)) / temp
    }
  }
  
  return(list(scores=gene.scores, diffs=mean.diff, cors=cors, wilcox=wilcox, corrected=wilcox.bh))
}

