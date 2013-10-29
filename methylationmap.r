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
	if (!is.matrix(tumors)) {
		tumors = matrixFromVector(tumors)
		rownames(tumors) = probes
	}
	if (!is.matrix(normals)) {
		normals = matrixFromVector(normals)
		rownames(normals) = probes
	}
	
	# Samples from tumors that have a matched normal
	tumors.matchedsamples = intersect(colnames(tumors), colnames(normals))
	# Rename the normal samples 
	colnames(normals) = paste(colnames(normals), "normal", sep="_")
	
	# Filter unknown probes
	selected.probes = intersect(intersect(rownames(tumors), rownames(normals)), probes)
	
	## Unpaired test
	tumors.groups = c(rep(1, times=ncol(tumors)), rep(2, times=ncol(normals)))
	names(tumors.groups) = c(colnames(tumors), colnames(normals))
	
	inpMat = cbind(tumors[selected.probes, , drop=FALSE], normals[selected.probes, , drop=FALSE])
	# Run un-paired Wilcoxon tests
	wilcox.p = doWilcox(inpMat=inpMat, groups=tumors.groups)
	
	## Paired test
	# Retain only the normal samples for which there is a tumor
	normals = normals[, tumors.matchedsamples, drop=FALSE]
		
	inpMat = cbind(tumors[selected.probes, tumors.matchedsamples, drop=FALSE], normals[selected.probes, , drop=FALSE])
	# Run paired Wilcoxon tests
	pairedwilcox.p = doWilcox(inpMat=inpMat, matchedSamples=tumors.matchedsamples)
	
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
# gene2probe named list mapping genes to probes. Names are gene symbols and values are vectors
# with probe accessions
# exprs.data should be a matrix with probes in rows and samples in columns
# meth.probes should be a character vector giving the probes to test
corMethExprs = function(meth.data, gene2probe, exprs.data, meth.probes) {
  # All genes contained in the expression data
  exprs.genes = rownames(exprs.data)
  commonSamples = intersect(colnames(meth.data), colnames(exprs.data))
	
  cors = numeric(length(meth.probes))
  n = character(length(meth.probes))
  j = 0
  for (gene in names(gene2probe)) {
	for (probe in intersect(gene2probe[[gene]], meth.probes)) {  
		j = j + 1
		cors[j] = cor(exprs.data[gene, commonSamples], 
					  meth.data[probe, commonSamples], 
					  method="spearman")
		n[j] = probe
    }
  }
  cors = cors[1:j]
  names(cors) = n[1:j]
  
  cors
}

# Calculates correlation between methylation probes and expression of the corresponding genes.
# This method deals with non-unique mappings of probes to genes by testing all
# of them.
# Correlation is calculated using all samples contained in both data matrices
# meth.data should be a matrix with probes in rows and samples in columns
# meth.ann has to be an annotation matrix mapping probes to genes;
# probe ids in 1st column, gene ids and gene regions with colnames "gene" and "gene_region", respectively.
# gene ids have to match ids used in the expression data.
# exprs.data should be a matrix with probes in rows and samples in columns
# meth.probes should be a character vector giving the probes to test
corMethExprsMult = function(meth.data, meth.ann, exprs.data, meth.probes, filt=c("all", "low", "high")) {
	filt = match.arg(filt)
	
	exprs.genes = rownames(exprs.data)
	samps = intersect(colnames(meth.data), colnames(exprs.data))
	
	# Correlations and their names
	cors = list()
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
		genes = unlist(strsplit(meth.ann[which(meth.ann[, 1] == x), "gene"], ";"))
		regions = unlist(strsplit(meth.ann[which(meth.ann[, 1] == x), "gene_region"], ";"))
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
	cors[[j]] = temp1[1:(i-1)]
	n[[j]] = temp2[1:(i-1)]
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

  # How many comparisons per gene exceed the threshold? 
  gene.scores = apply(methylation, 1, function(x) { length(which(x > threshold)) })	
  # Is there at least one comparison exceeding the threshold?
  if (score)
    gene.scores = sapply(gene.scores, function(x) { as.integer(x > 0) }) 
  
  return(gene.scores)
}

##' Preprocessing of methylation annotation. Select all probes that map to
##' any of the genes of interest and map these genes to their probes.
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param genes vector with gene symbols 
##' @param probe.annotation methylation probe annotation matrix. The first column has to be the 
##' probe id, "gene" column giving gene IDs (possibily separated by ";")
##' @return named list with two entries ("selected.probes" and "gene2probe")
##' @author Andreas Schlicker
genesAndProbes = function(tumors, normals, genes, probe.annotation) {
	require(stringr) || stop("Could not load package \"stringr\"!")
	
	common.probes = intersect(rownames(tumors), rownames(normals))
	
	# Does the probe in the given row map to any of the genes of interest?
	found.gene = unlist(lapply(lapply(probe.annotation[common.probes, "gene"], 
									  function(x) { unlist(unique(str_split(x, ";"))) } ), 
							   function(y) { any(y %in% genes) }))
	
	# All probes that map to any of the genes
	selected.probes = intersect(common.probes, 
			rownames(probe.annotation)[found.gene])
	
	# Mapping of genes to the probes
	gene2probe = list()
	for (probe in selected.probes) {
		temp = unique(str_split(probe.annotation[probe, "gene"], ";")[[1]])
		for (gene in temp) { 
			if (is.null(gene2probe[[gene]])) {
				gene2probe[[gene]] = c()
			}
			gene2probe[[gene]] = c(gene2probe[[gene]], probe)
		}
	}
	
	list(selected.probes=selected.probes, gene2probe=gene2probe)
}

##' Performs all steps of the methylation analysis. First, probes that don't meet the
##' "diff.cutoff" are filtered out. Second, probes without significant correlation
##' with gene expression are filtered out. Third, probes without significant difference
##' in methylation are filtered out. 
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param genes vector with gene symbols 
##' @param exprs gene expression matrix with genes in rows and samples in columns
##' @param probe.annotation methylation probe annotation matrix. The first column has to be the 
##' probe id, "gene" column giving gene IDs (possibily separated by ";"), "gene_region" column
##' giving name of gene region for every gene (possibily separated by ";"), "chrom" and
##' "pos" giving the position.
##' @param selected.probes vector with all probes to investigate (as returned by genesAndProbes())
##' @param gene2probe mapping of genes to their methylation probes (as returned by genesAndProbes())
##' @param samples vector with samples to include in the analysis, if == NULL, all samples are used; default: NULL
##' @param wilcox.cutoff significance cut-off for Wilcoxon test; this cut-off is applied to Benjamini-Hochberg-
##' corrected p-values; default: 0.05
##' @param diff.cutoff minimum methylation difference for a probe to be considered; default: 0.1
##' @param regulation either "down" or "up", indicating whether down- or upregulation of the genes are tested.
##' Downregulation implies a gain of methylation in promoter regions. default: down
##' @param cor.min minimum correlation that a probe must exhibit to be considered significant; default: -1.0
##' @param cor.max maximum correlation that a probe must exhibit to be considered significant; default: -0.1
##' @param gene.region boolean indicating whether gene region information should be considered; default: FALSE.
##' If this is FALSE, correlation for all probes is tested to be in the interval [cor.min; cor.max]. If this is
##' TRUE, probes in the gene body are treated in the opposite way as all other probes because they are expected
##' to have a positive correlation with gene expression.
##' @param stddev used to assess whether a samples is affected by a change in methylation; the sample has to deviate
##' more than this many standard deviations in normal samples from the mean in normal samples; default: 1  
##' @return named list with ratio of significant probes for each genes, difference in average methylation between tumors and
##' normals, correlation values between methylation and expression, paired and unpaired Wilcoxon p-values and corrected
##' Wilcoxon p-values
##' @author Andreas Schlicker
doMethylationAnalysis = function(tumors, 
								 normals, 
								 genes, 
								 exprs, 
								 probe.annotation,
								 selected.probes,
								 gene2probe,
								 samples=NULL, 
								 wilcox.cutoff=0.05, 
								 diff.cutoff=0.1, 
								 regulation=c("down", "up"), 
								 cor.min=-1.0, 
								 cor.max=-0.1,
								 gene.region=FALSE,
								 stddev=1) {
							
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
  
  # Mean methylation filtering
  mean.tumor = apply(tumors[selected.probes, tumor.samples, drop=FALSE], 1, mean, na.rm=TRUE)
  mean.normal = apply(normals[selected.probes, normal.samples, drop=FALSE], 1, mean, na.rm=TRUE)
  mean.diff = mean.tumor - mean.normal
  selected.probes = names(mean.diff[compare(mean.diff, diff.cutoff)])
  selected.probes = selected.probes[!is.na(selected.probes)]
  
  # Correlation between methylation and expression data filtering
  cors = c()
  if (length(selected.probes) > 0) {
    temp = tumors[selected.probes, tumor.samples, drop=FALSE]
    # gene.region determines whether gene region information is available for the
	# methylation probes. If so, body probes are treated differently from other probes.
	if (!gene.region) {
    	cors = corMethExprs(temp, gene2probe, exprs, selected.probes)
    	selected.probes = names(cors[which(cors <= cor.max & cors >= cor.min)])
    } else {
		cors = corMethExprsMult(temp, probe.annotation, exprs, selected.probes)
		# The correlation names have the form <meth.probe_gene_gene.region>
		# Split everything up into matrix form:
		# meth.probe, gene, gene.region, correlation
		cors = cbind(data.frame(do.call("rbind", strsplit(names(cors), "_"))), cors)
		colnames(cors) = c("meth.probe", "gene", "gene.region", "cor.val")
		cors[, "meth.probe"] = as.character(cors[, "meth.probe"])
		selected.probes = union(cors[which(cors[, "gene.region"] != "Body" & cors[, "cor.val"] <= cor.max & cors[, "cor.val"] >= cor.min), "meth.probe"],
								cors[which(cors[, "gene.region"] == "Body" & cors[, "cor.val"] <= (-1*cor.min) & cors[, "cor.val"] >= (-1*cor.max)), "meth.probe"])
	}
  }
  
  wilcox = list(unpaired=c(), paired=c())
  wilcox.bh = list(unpaired=c(), paired=c())
  selected.probes.unpaired = c()
  selected.probes.paired = c()
  affected.samples.unpaired = list()
  affected.samples.paired = list()
  if (length(selected.probes) > 0) {
    tt = tumors[selected.probes, tumor.samples, drop=FALSE]
    tn = normals[selected.probes, normal.samples, drop=FALSE]
    
    wilcox = runMethComp(tt, tn, selected.probes)
    wilcox.bh = list(unpaired=p.adjust(wilcox$unpaired, method="BH"), paired=p.adjust(wilcox$paired, method="BH"))
    selected.probes.unpaired = names(wilcox.bh$unpaired[which(wilcox.bh$unpaired < wilcox.cutoff)])
    selected.probes.paired = names(wilcox.bh$paired[which(wilcox.bh$paired < wilcox.cutoff)])
	# Check which samples are affected by the difference in methylation. 
	nonbody = countAffectedSamples(intersect(selected.probes.unpaired, cors[which(cors[, "gene.region"] != "Body"), "meth.probe"]), 
													 tumors[selected.probes.unpaired, , drop=FALSE], 
													 normals[selected.probes.unpaired, , drop=FALSE], 
													 regulation=switch(regulation, up="down", down="up"), 
													 stddev, FALSE)
	body = countAffectedSamples(intersect(selected.probes.unpaired, cors[which(cors[, "gene.region"] == "Body"), "meth.probe"]), 
													 tumors[selected.probes.unpaired, , drop=FALSE], 
													 normals[selected.probes.unpaired, , drop=FALSE], 
													 regulation=regulation, 
													 stddev, FALSE)
	affected.samples.unpaired = list(summary=rbind(nonbody$summary, body$summary), samples=c(nonbody$samples, body$samples))
	nonbody = countAffectedSamples(intersect(selected.probes.paired, cors[which(cors[, "gene.region"] != "Body"), "meth.probe"]), 
													 tumors[selected.probes.paired, , drop=FALSE], 
													 normals[selected.probes.paired, , drop=FALSE], 
						 							 regulation=switch(regulation, up="down", down="up"), 
													 stddev, TRUE)
	body = countAffectedSamples(intersect(selected.probes.paired, cors[which(cors[, "gene.region"] == "Body"), "meth.probe"]), 
													 tumors[selected.probes.paired, , drop=FALSE], 
													 normals[selected.probes.paired, , drop=FALSE], 
													 regulation=regulation, 
													 stddev, TRUE)
	affected.samples.paired = list(summary=rbind(nonbody$summary, body$summary), samples=c(nonbody$samples, body$samples))
  }
  
  gene.scores = list(unpaired=rep(0.0, length(genes)), paired=rep(0.0, length(genes))) 
  names(gene.scores$unpaired) = genes
  names(gene.scores$paired) = genes
  samples = list(unpaired=list(), paired=list())
  summary = list(unpaired=matrix(0, ncol=2, nrow=length(genes)), paired=matrix(0, ncol=2, nrow=length(genes)))
  rownames(summary$unpaired) = genes
  colnames(summary$unpaired) = c("absolute", "relative")
  rownames(summary$paired) = genes
  colnames(summary$paired) = c("absolute", "relative")
  nTumors = ncol(tumors)
  for (gene in genes) {
    # Get all probes for that gene
    geneProbes = gene2probe[[gene]]
	allProbes = rownames(probe.annotation[geneProbes, , drop=FALSE])
    
    temp = length(geneProbes)
    if (temp > 0) {
      gene.scores$unpaired[gene] = length(intersect(allProbes, selected.probes.unpaired)) / temp
      gene.scores$paired[gene] = length(intersect(allProbes, selected.probes.paired)) / temp
	  samples$unpaired[[gene]] = unique(unlist(affected.samples.unpaired$samples[intersect(allProbes, selected.probes.unpaired)]))
	  samples$paired[[gene]] = unique(unlist(affected.samples.paired$samples[intersect(allProbes, selected.probes.paired)]))
	  summary$unpaired[gene, "absolute"] = length(samples$unpaired[[gene]])
	  summary$unpaired[gene, "relative"] = summary$unpaired[gene, "absolute"] / nTumors
	  summary$paired[gene, "absolute"] = length(samples$paired[[gene]])
	  summary$paired[gene, "relative"] = summary$paired[gene, "absolute"] / nTumors
    }
  }
  
  return(list(scores=gene.scores, diffs=mean.diff, cors=cors, wilcox=wilcox, corrected=wilcox.bh, samples=samples, summary=summary))
}

