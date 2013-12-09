##' Preprocessing of methylation annotation. Select all probes that map to
##' any of the genes of interest and map these genes to their probes.
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param genes vector with gene symbols; default: NULL (all genes)
##' @param probe.annotation methylation probe annotation matrix. The first column has to be the 
##' probe id, "gene" column giving gene IDs (possibily separated by ";")
##' @param gene2probe nested named list that maps probes to regions of genes; names of the outer list
##' map to gene IDs used in "genes"; names of inner lists are gene regions
##' @param regions vector with gene regions to remove; default: "" to keep probes mapping to any region
##' @param snps probes with SNPs in any of these locations will be removed; default: any SNP in the probe
##' or target region
##' @return named list with two entries ("selected.probes" and "gene2probe")
##' @author Andreas Schlicker
filterProbes = function(tumors, normals, 
						genes=NULL, probe.annotation, gene2probe, probe2gene,
						regions=c(""), 
						snps=c("SNP_target", "SNP_within_10", "SNP_outside_10", "SNP_probe")) {
	require(stringr) || stop("Could not load required package \"stringr\"!")
	
	# Probes in the two data set that also have annotation
	common.probes = intersect(rownames(tumors), intersect(rownames(normals), rownames(probe.annotation)))
	
	# If genes == NULL, all genes
	if (is.null(genes)) {
		genes = names(gene2probe)
	}
	# Get probe to gene and region mapping for all genes
	res.g2p = gene2probe[genes]
	
	## Filter probes by gene region
	# If analysis should be restricted to some gene regions, remove all others 
	if (length(regions) > 1) {
		for (x in 1:length(res.g2p)) {
			res.g2p[[x]][regions] = NULL
		}
	}
	# Vector with probes to remove
	excl.probes = setdiff(unique(unlist(res.g2p)), common.probes)
	
	## Filter probes by SNPS
	if (!is.null(snps) && length(snps) > 0) {
		# We need to filter out probes containing SNPs in the given categories
		# Find all probes that have a SNP in at least one category
		excl.snps = apply(infinium450.probe.ann[unique(unlist(res.g2p)), snps, drop=FALSE], 1, function(x) { any(x) })
		excl.snps = names(excl.snps)[excl.snps]
		excl.probes = union(excl.probes, excl.snps)
	}
	
	if (length(excl.probes) > 0) {
		res.selprobes = setdiff(unique(unlist(res.g2p)), excl.probes)
		# Clean up the gene to probe mapping as well
		for (p in excl.probes[1:10000]) {
			for (region in names(probe2gene[[p]])) {
				for (gene in probe2gene[[p]][[region]]) {
					temp = res.g2p[[gene]][[region]]
					temp = temp[-match(p, temp)]
					if (length(temp) == 0) {
						res.g2p[[gene]][[region]] = NULL
					} else {
						res.g2p[[gene]][[region]] = temp
					}
				}
			}
		}
	}
	
	res.g2p.flat = vector("list", length(res.g2p))
	names(res.g2p.flat) = names(res.g2p)
	for (i in 1:length(res.g2p)) {
		res.g2p.flat[[i]] = unlist(res.g2p[[i]])
	}
	
	list(selected.probes=res.selprobes, gene2probe=res.g2p, gene2probe.flat=res.g2p.flat)
}

##' Calculates Spearman correlation between methylation probes and expression of the corresponding genes and 
##' the corresponding p-values. Correlation is calculated using all samples contained in both data matrices.
##' @param meth.data matrix with probes in rows and samples in columns. Correlations will be computed 
##' for all probes contained in this matrix (identified by rownames).
##' @param probe2gene.flat named list mapping probes to genes
##' @param exprs.data matrix with probes in rows and samples in columns
##' @return data.frame with 4 columns: probe, gene, correlation, p.value
##' @author Andreas Schlicker
corMethExprs = function(meth.data, probe2gene.flat, exprs.data) {
	# Samples with methylation and expression data
	samps = intersect(colnames(meth.data), colnames(exprs.data))
	exprs.data = exprs.data[, samps, drop=FALSE]
	meth.data = meth.data[, samps, drop=FALSE]
	
	# Correlations
	res = list()
	# Length of temporary vectors
	LENGTH = 1000
	# Temporary vectors
	temp1 = data.frame(probe=character(LENGTH), gene=character(LENGTH), 
					   cor=rep(NA, times=LENGTH), p.value=rep(NA, times=LENGTH), 
					   stringsAsFactors=FALSE)
	# Counter for current element in temporary vector
	i = 1
	k = 1
	testProbes = rownames(meth.data)
	for (j in 1:length(testProbes)) {
		tmpMeth = meth.data[j, ]
		x = testProbes[j]
		
		# Get all genes this probe maps to
		genes = probe2gene.flat[[x]]
		
		# No space left in temporary data.frame
		if ((i+length(genes)) > LENGTH) {
			res[[k]] = temp1[1:(i-1), ]
			k = k + 1
			i = 1
		}
	
		for (gene in genes) {
			tempCor = tryCatch(cor.test(exprs.data[gene, ], 
							   			meth.data[x, ], 
							   			method="spearman", 
							   			use="pairwise.complete.obs",
							   			exact=FALSE),
							   error=function(e) NA)
			if (class(tempCor) == "htest") {
				temp1[i, ] = c(x, gene, tempCor$estimate, tempCor$p.value)
				i = i + 1
			}
		}
	}
	# Save the last values
	res[[k]] = temp1[1:(i-1), ]
	res = do.call("rbind", res)
	colnames(res) = c("probe", "gene", "cor", "cor.p")
	res[, "cor"] = as.numeric(res[, "cor"])
	res[, "cor.p"] = as.numeric(res[, "cor.p"])
	
	res
}

##' Performs all steps of the methylation analysis. First, probes that don't meet the
##' "diff.cutoff" are filtered out. Second, probes without significant correlation
##' with gene expression are filtered out. Third, probes without significant difference
##' in methylation are filtered out. 
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param exprs gene expression matrix with genes in rows and samples in columns
##' @param probe2gene.flat named list mapping probes to genes 
##' @param selected.probes vector with all probes to investigate (as returned by genesAndProbes()); 
##' default: NULL (test all probes in both tumors and normals) 
##' @param samples vector with samples to include in the analysis, if == NULL, all samples are used; default: NULL
##' @return named list with three entries: diff: difference in average methylation between tumors and
##' normals; cors: correlation values between methylation and expression, with p-values and FDR;
##' wilcox: paired and unpaired Wilcoxon p-values and FDRs for each probe
##' @param paired boolean indicating whether doing paired or unpaired analysis; default: FALSE
##' @author Andreas Schlicker
doMethylationAnalysis = function(tumors, 
								 normals, 
								 exprs, 
								 probe2gene.flat,
								 selected.probes=NULL,
								 samples=NULL, 
								 paired=FALSE) {
	common.probes = doFilter(rownames(tumors), rownames(normals), selected.probes, TRUE)
	if (length(common.probes[[1]]) == 0) {
		stop("doMethylationAnalysis: No probe of interest is contained in methylation data of both tumors and normals!")
	}
	
	filtered.samples = doFilter(colnames(tumors), colnames(normals), samples, paired)
	
	if (paired && (length(filtered.samples[[1]]) == 0 || length(filtered.samples[[2]]) == 0)) {
		paired = FALSE
		filtered.samples = doFilter(colnames(tumors), colnames(normals), samples, FALSE)
		warning("doMethylationAnalysis: No paired expression samples found. Performing unpaired analysis!")
	}
	
	tumors = tumors[common.probes[[1]], filtered.samples[[1]], drop=FALSE]
  	normals = normals[common.probes[[2]], filtered.samples[[2]], drop=FALSE]
  
  	# Mean methylation filtering
  	mean.diff = meanDiff(tumors, normals)
  
  	# Correlation between methylation and expression data
	cors = corMethExprs(tumors, probe2gene.flat, exprs)
	
	# Wilcoxon test
	wilcox = doWilcox(tumors, normals, paired)
	
	list(diffs=mean.diff, cors=cors, wilcox=wilcox)
}


##' Calculates the ratio of significant probes per gene and the list of samples that are affected by changes in these probes.
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param meth.analysis list returned by doMethylationAnalysis()
##' @param gene2probe.flat named list mapping probes to genes and gene regions
##' @param probe2gene nested list mapping genes to probes and gene regions
##' @param genes character vector of genes; default: NULL (test all genes in tumors)
##' @param wilcox.FDR significance cut-off for Wilcoxon FDR; default=0.05
##' @param cor.FDR significance cut-off for correlation FDR; default=0.05
##' @param diff.cutoff methylation difference needs to be smaller or greater than this cut-off to be considered significant; default=0.1
##' @param regulation either "down" or "up" for finding genes that are regulated in the corresponding direction; default="down"
##' @param gene.region boolean indicating whether gene region annotation should be taken into account; default=TRUE
##' @param stddev how many standard deviations does a sample have to be away from the mean to be considered affected; default=1
##' @param score.cutoff double, a gene scores 1, if it exceeds this cut-off
##' @return named list with scores for genes ("scores"), number of affected samples ("summary") and the lists of affected samples ("samples")
##' @author Andreas Schlicker
summarizeMethylation = function(tumors,
				 				normals,
							 	meth.analysis,
							 	gene2probe.flat,
							 	probe2gene,
							 	genes=NULL,
							 	wilcox.FDR=0.05, 
							 	cor.FDR=0.05,
							 	diff.cutoff=0.1, 
							 	regulation=c("down", "up"),
							 	gene.region=TRUE,
							 	stddev=1,
								score.cutoff=0) { 
	regulation = match.arg(regulation)
	
	# Get the correct comparison function
	# If we want to find genes with higher methylation in tumors get the greaterThan function
	# If we want to find genes with lower methylation in tumors, get the smallerThan function
	compare = switch(regulation, down=greaterThan, up=smallerThan)
	# Gene body probes should behave differently
	compareBody = switch(regulation, down=smallerThan, up=greaterThan)
	
	# Restrict the analysis to probes that match to the genes of interest
	if (!is.null(genes)) {
		significant.probes = doFilter(rownames(tumors), rownames(normals), unlist(gene2probe.flat[genes]), TRUE)[[1]]
	} else {
		significant.probes = doFilter(rownames(tumors), rownames(normals), NULL, TRUE)[[1]]
		genes = names(gene2probe.flat)
	}
	all.probes = significant.probes
	
	meth.analysis$cors = meth.analysis$cors[which(meth.analysis$cors[, "probe"] %in% significant.probes), ]
	meth.analysis$diffs = meth.analysis$diffs[which(names(meth.analysis$diffs) %in% significant.probes)]
	meth.analysis$wilcox = meth.analysis$wilcox[which(names(meth.analysis$wilcox) %in% significant.probes)]
	
	# Correlation significance filter
	# Multiple testing correction
	meth.analysis$cors = cbind(meth.analysis$cors, cor.FDR=p.adjust(meth.analysis$cors[, "cor.p"], method="BH"))
	significant.cors = meth.analysis$cors[which(meth.analysis$cors[, "cor.FDR"] <= cor.FDR), ]
	
	tsc.body = integer(nrow(significant.cors))
	bodyCounter = 1
	tsc.rest = integer(nrow(significant.cors))
	restCounter = 1
	if (!gene.region) {
		tsc.rest = which(significant.cors[, "cor"] < 0)
	} else {
		tmpSigCorsGreater = significant.cors[, "cor"] > 0
		tmpSigCorsSamller = significant.cors[, "cor"] < 0
		for (i in 1:nrow(significant.cors)) {
			if (significant.cors[i, "gene"] %in% probe2gene[[significant.cors[i, "probe"]]][["Body"]]) {
				# We have a body probe --> positive correlation is expected
				# Note: This might not be true for 1st exon probes, but they are annotated separately
				if (tmpSigCorsGreater[i]) {
					tsc.body[bodyCounter] = i
					bodyCounter = bodyCounter + 1
				}
			} else if (tmpSigCorsSamller[i]) {
				# We have a probe that is not in the gene body --> expect negative correlation
				tsc.rest[restCounter] = i
				restCounter = restCounter + 1
			}
		}
	}
	significant.probes = intersect(significant.probes, unique(significant.cors[c(tsc.body[1:(bodyCounter-1)], tsc.rest[1:(restCounter-1)]), "probe"]))
	
	# Calculate multiple testing correction for remaining probes
	meth.analysis$wilcox = cbind(wilcox.p=meth.analysis$wilcox, wilcox.FDR=NA)
	meth.analysis$wilcox[significant.probes, "wilcox.FDR"] = p.adjust(meth.analysis$wilcox[significant.probes, "wilcox.p"], method="BH")
	# Mean difference filter 
	significant.probes = intersect(significant.probes, 
								   union(names(which(compare(meth.analysis$diffs[unique(significant.cors[tsc.rest, "probe"])], diff.cutoff))),
										 names(which(compareBody(meth.analysis$diffs[unique(significant.cors[tsc.body, "probe"])], diff.cutoff)))))
	# Wilcoxon significance filter
	significant.probes = names(which(meth.analysis$wilcox[significant.probes, "wilcox.FDR"] <= wilcox.FDR))
	
	## Find affected samples per probe
	bodyprobes = intersect(significant.probes, significant.cors[tsc.body, "probe"])
	nonbodyprobes = intersect(significant.probes, significant.cors[tsc.rest, "probe"])
	nonbody = countAffectedSamples(all.probes,
								   nonbodyprobes, 
								   tumors[nonbodyprobes, , drop=FALSE], 
								   normals[nonbodyprobes, , drop=FALSE], 
								   regulation=switch(regulation, up="down", down="up"), 
								   stddev, 
								   paired=FALSE)
	body = countAffectedSamples(all.probes,
								bodyprobes,
								tumors[bodyprobes, , drop=FALSE], 
								normals[bodyprobes, , drop=FALSE], 
								regulation=regulation, 
								stddev, 
								paired=FALSE)
	
	## Summarize results
	nTumors = ncol(tumors)
	gene.scores = double(length(genes))
	names(gene.scores) = genes
	summary = matrix(0, ncol=2, nrow=length(genes))
	rownames(summary) = genes
	colnames(summary) = c("absolute", "relative")
	samples = list()
	for (gene in genes) {
		allProbes = gene2probe.flat[[gene]]
		# Get the ratio of (#significant probes with higher methylation) / (#all probes)
		gene.scores[gene] = length(intersect(significant.probes, allProbes)) / length(allProbes)
		
		samples[[gene]] = union(unlist(body$samples[intersect(allProbes, names(body$samples))]),
								unlist(nonbody$samples[intersect(allProbes, names(nonbody$samples))]))
		summary[gene, "absolute"] = length(samples[[gene]])
		summary[gene, "relative"] = summary[gene, "absolute"] / nTumors
	}
	
	gene.scores = sapply(gene.scores, function(x) { as.integer(x > score.cutoff) })
	
	list(scores=gene.scores, summary=summary, samples=samples, meth.analysis=meth.analysis)
}

##' Summarizes the results of analyzing methylation data.
##' @param methylation matrix containing the results of at least one methylation set. 
##' Typically this will be the output of the "summarizeMethExprs" function.
##' A gene receives a score of one if at least on entry for that gene exceeds the
##' given threshold.
##' @param genes a character vector containing the IDs of the genes to analyze; IDs have
##' to map to row names in the "methylation" matrix
##' @param threshold floating point number giving the threshold that has to be met;
##' defaults to 0
##' @param score either "atleast" (default) or "absolute". If set to "atleast", a gene
##' receives a score of 1 if at least one column exceeds the threshold. If set to
##' "absolute", the score of the gene will be equal to the number of columns 
##' exceeding the threshold.
##' @return numeric vector with the gene scores
##' @author Andreas Schlicker
combineMethylation = function(methylation, genes, threshold=0, score=c("atleast", "absolute")) {
	score = match.arg(score)
	score = switch(score, atleast=TRUE, absolute=FALSE)
	
	# How many comparisons per gene exceed the threshold? 
	gene.scores = apply(methylation, 1, function(x) { length(which(x > threshold)) })	
	# Is there at least one comparison exceeding the threshold?
	if (score) {
		gene.scores = sapply(gene.scores, function(x) { as.integer(x > 0) })
	}
	
	gene.scores
}
