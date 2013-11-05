##' Preprocessing of methylation annotation. Select all probes that map to
##' any of the genes of interest and map these genes to their probes.
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param genes vector with gene symbols
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
		genes, probe.annotation, gene2probe,
		regions=c(""), 
		snps=c("SNP_target", "SNP_within_10", "SNP_outside_10", "SNP_probe")) {
	require(stringr) || stop("Could not load required package \"stringr\"!")
	
	# Probes in the two data set that also have annotation
	common.probes = intersect(rownames(tumors), intersect(rownames(normals), rownames(probe.annotation)))
	
	# Get probe to gene and region mapping for all genes 
	res.g2p = gene2probe[genes]
	# If analysis should be restricted to some gene regions, remove all others 
	if (regions != "") {
		for (x in names(res.g2p)) {
			for (n in intersect(names(res.g2p[[x]]), regions)) {
				res.g2p[[x]][[n]] = NULL 
			}
		}
	}
	
	# All probes that we look at
	res.selprobes = intersect(unique(unlist(unlist(res.g2p))), common.probes)
	
	# We need to filter out probes containing SNPs in the given categories
	if (length(snps) > 0) {
		# Find all probes that have a SNP in at least one category
		excl.snps = apply(infinium450.probe.ann[res.selprobes, snps, drop=FALSE], 1, function(x) { any(x) })
		excl.snps = names(excl.snps)[excl.snps]
		
		# That's the remining probes
		res.selprobes = setdiff(res.selprobes, excl.snps)
		# Clean up the gene to probe mapping as well
		for (p in excl.snps) {
			# All genes and their regions with the current probe
			gene_regions = str_split(probe.annotation[p, "Symbol_region_unique"], ";")[[1]]
			for (gene_region in gene_regions) {
				# Get the gene and the region, remove the probe; if there is no probe left, delete the region
				gs = str_split(gene_region, "_")[[1]]
				gene = gs[1]
				region = gs[2]
				
				res.g2p[[gene]][[region]] = setdiff(res.g2p[[gene]][[region]], p)
				if (length(res.g2p[[gene]][[region]]) == 0) {
					res.g2p[[gene]][[region]] == NULL
				}
			}
		}
	}
	
	list(selected.probes=res.selprobes, gene2probe=res.g2p)
}

##' Compares methylation data of tumor and normal samples using unpaired Wilcoxon test.
##' @param tumors methylation matrix with samples in the columns and probes in the rows
##' @param normals methylation matrix with samples in the columns and probes in the rows
##' @return data.frame with two columns containing results from paired and un-paired tests
##' @author Andreas Schlicker
runMethComp = function(tumors, normals) {
	# Samples from tumors that have a matched normal
	tumors.matchedsamples = intersect(colnames(tumors), colnames(normals))
	# Rename the normal samples 
	colnames(normals) = paste(colnames(normals), "normal", sep="_")
	
	## Unpaired test
	tumors.groups = c(rep(1, times=ncol(tumors)), rep(2, times=ncol(normals)))
	names(tumors.groups) = c(colnames(tumors), colnames(normals))
	# Run un-paired Wilcoxon tests
	doWilcox(inpMat=cbind(tumors, normals), groups=tumors.groups)
}

##' Calculates Spearman correlation between methylation probes and expression of the corresponding genes and 
##' the corresponding p-values. Correlation is calculated using all samples contained in both data matrices.
##' @param meth.data matrix with probes in rows and samples in columns. Correlations will be computed 
##' for all probes contained in this matrix (identified by rownames).
##' @param meth.ann annotation matrix mapping probes to genes; rownames of the matrix need to be probe IDs,
##' and column "Gene_symbols_unique" needs to contain all genes it maps to; gene IDs have to match IDs used 
##' in the expression data.
##' @param exprs.data matrix with probes in rows and samples in columns
##' @return data.frame with 4 columns: probe, gene, correlation, p.value
##' @author Andreas Schlicker
corMethExprs = function(meth.data, meth.ann, exprs.data) {
	# Genes with expression data
	exprs.genes = rownames(exprs.data)
	# Samples with methylation and expression data
	samps = intersect(colnames(meth.data), colnames(exprs.data))
	
	# Correlations
	res = list()
	# Length of temporary vectors
	LENGTH = 100
	# Temporary vectors
	temp1 = data.frame(probe=character(LENGTH), gene=character(LENGTH), 
					   cor=double(LENGTH), p.value=double(LENGTH), 
					   stringsAsFactors=FALSE)
	# Counter for current element in temporary vector
	i = 1
	
	for (x in rownames(meth.data)) {
		# No space left in temporary list
		if (i > LENGTH) {
			cors = c(cors, temp1)
			temp1 = data.frame(probe=character(LENGTH), gene=character(LENGTH), 
							   cor=double(LENGTH), p.value=double(LENGTH), 
							   stringsAsFactors=FALSE)
			i = 1
		}
		
		# Get all genes and gene regions and combine them
		genes = intersect(unlist(strsplit(meth.ann[x, "Gene_symbols_unique"], ";")),
						  exprs.genes)
		for (gene in genes) {
			tempCor = cor.test(exprs.data[gene, samps], meth.data[x, samps], method="spearman", exact=FALSE)
			temp1[i, ] = c(x, gene, tempCor$estimate, tempCor$p.value)
		}
	}
	# Save the last values
	cors = c(cors, temp1[1:(i-1), ])
	
	res = do.call("rbind", cors)
	colnames(res) = c("probe", "gene", "cor", "cor.p")
}

##' Performs all steps of the methylation analysis. First, probes that don't meet the
##' "diff.cutoff" are filtered out. Second, probes without significant correlation
##' with gene expression are filtered out. Third, probes without significant difference
##' in methylation are filtered out. 
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param exprs gene expression matrix with genes in rows and samples in columns
##' @param probe.annotation methylation probe annotation matrix. 
##' @param selected.probes vector with all probes to investigate (as returned by genesAndProbes())
##' @param samples vector with samples to include in the analysis, if == NULL, all samples are used; default: NULL
##' @return named list with three entries: diff: difference in average methylation between tumors and
##' normals; cors: correlation values between methylation and expression, with p-values and FDR;
##' wilcox: paired and unpaired Wilcoxon p-values and FDRs for each probe
##' @author Andreas Schlicker
doMethylationAnalysis = function(tumors, 
								 normals, 
								 exprs, 
								 probe.annotation,
								 selected.probes,
								 samples=NULL) {
	if (!is.null(samples)) {
	  	tumors = tumors[, intersect(samples, colnames(tumors)), drop=FALSE]
    	normals = normals[, intersect(samples, colnames(normals)), drop=FALSE]
  	}
  
  	# Restrict to selected probes
	common.probes = intersect(selected.probes, intersect(rownames(tumors), rownames(normals)))
	if (length(common.probes) == 0) {
		stop("No probes to test in doMethylationAnalysis!")
	}
  	tumors = tumors[common.probes, , drop=FALSE]
  	normals = normals[common.probes, , drop=FALSE]
  
  	# Mean methylation filtering
  	mean.diff = meanDiff(tumors, normals)
  
  	# Correlation between methylation and expression data
  	cors = corMethExprs(tumors, probe.annotation, exprs)
	
	# Wilcoxon test
	wilcox = runMethComp(tumors, normals)
	
	list(diffs=mean.diff, cors=cors, wilcox=wilcox)
	
	# gene.region determines whether gene region information is available for the
	# methylation probes. If so, body probes are treated differently from other probes.
	#if (!gene.region) {
	#	cors = corTestMethExprs(temp, gene2probe, exprs, selected.probes)
	#	pvalues = unlist(lapply(cors, function(x) { x$p.value }))
	#	fdrs = p.adjust(pvalues, method="BH")
	#	cors = cbind(data.frame(do.call("rbind", strsplit(names(cors), "_"))), 
	#				 "any",
	#				 unlist(lapply(cors, function(x) { x$estimate })),
	#				 pvalues, fdrs)
	#   colnames(cors) = c("meth.probe", "gene", "gene.region", "cor.val", "pvalue", "FDR")
	#	rownames(cors) = str_replace(rownames(cors), "\\.rho", "")
	#		 
	#	#selected.probes = names(cors[which(cors <= cor.max & cors >= cor.min)])
	#	selected.probes = unique(names(fdrs)[fdrs < wilcox.cutoff])
    #} else {
	#	cors = corTestMethExprsMult(temp, probe.annotation, exprs, selected.probes)
	#	# The correlation names have the form <meth.probe_gene_gene.region>
	#	# Split everything up into matrix form:
	#	# meth.probe, gene, gene.region, correlation
	#	pvalues = unlist(lapply(cors, function(x) { x$p.value }))
	#	fdrs = p.adjust(pvalues, method="BH")
	#	cors = cbind(data.frame(do.call("rbind", strsplit(names(cors), "_"))), 
	#				 unlist(lapply(cors, function(x) { x$estimate })),
	#				 pvalues, fdrs)
	#	colnames(cors) = c("meth.probe", "gene", "gene.region", "cor.val", "pvalue", "FDR")
	#	rownames(cors) = str_replace(rownames(cors), "\\.rho", "")
	#	cors[, "meth.probe"] = as.character(cors[, "meth.probe"])
	#	#selected.probes = union(cors[which(cors[, "gene.region"] != "Body" & cors[, "cor.val"] <= cor.max & cors[, "cor.val"] >= cor.min), "meth.probe"],
	#	#						cors[which(cors[, "gene.region"] == "Body" & cors[, "cor.val"] <= (-1*cor.min) & cors[, "cor.val"] >= (-1*cor.max)), "meth.probe"])
	#	selected.probes = union(cors[which(cors[, "gene.region"] != "Body" & cors[, "cor.val"] > 0 & cors[, "pvalue"] < wilcox.cutoff), "meth.probe"],
	#							cors[which(cors[, "gene.region"] == "Body" & cors[, "cor.val"] < 0 & cors[, "pvalue"] < wilcox.cutoff), "meth.probe"])
	#	#selected.probes = cors[which(cors[, "FDR"] < wilcox.cutoff), "meth.probe"]
	#}
  	#}
  
#	if (length(selected.probes) > 0) {
#	    tt = tumors[selected.probes, tumor.samples, drop=FALSE]
#    	tn = normals[selected.probes, normal.samples, drop=FALSE]
#    
#	    wilcox = runMethComp(tt, tn, selected.probes)
#    	wilcox.bh = list(unpaired=p.adjust(wilcox$unpaired, method="BH"), paired=p.adjust(wilcox$paired, method="BH"))
#    	selected.probes.unpaired = names(wilcox.bh$unpaired[which(wilcox.bh$unpaired < wilcox.cutoff)])
#    	selected.probes.paired = names(wilcox.bh$paired[which(wilcox.bh$paired < wilcox.cutoff)])
#		# Check which samples are affected by the difference in methylation. 
#		nonbody = countAffectedSamples(intersect(selected.probes.unpaired, cors[which(cors[, "gene.region"] != "Body"), "meth.probe"]), 
#												 tumors[selected.probes.unpaired, , drop=FALSE], 
#											     normals[selected.probes.unpaired, , drop=FALSE], 
#											 	 regulation=switch(regulation, up="down", down="up"), 
#												 stddev, FALSE)
#		body = countAffectedSamples(intersect(selected.probes.unpaired, cors[which(cors[, "gene.region"] == "Body"), "meth.probe"]), 
#			  								  tumors[selected.probes.unpaired, , drop=FALSE], 
#											  normals[selected.probes.unpaired, , drop=FALSE], 
#											  regulation=regulation, 
#										  	  stddev, FALSE)
#		affected.samples.unpaired = list(summary=rbind(nonbody$summary, body$summary), samples=c(nonbody$samples, body$samples))
#		nonbody = countAffectedSamples(intersect(selected.probes.paired, cors[which(cors[, "gene.region"] != "Body"), "meth.probe"]), 
#													 tumors[selected.probes.paired, , drop=FALSE], 
#													 normals[selected.probes.paired, , drop=FALSE], 
#						 							 regulation=switch(regulation, up="down", down="up"), 
#													 stddev, TRUE)
#		body = countAffectedSamples(intersect(selected.probes.paired, cors[which(cors[, "gene.region"] == "Body"), "meth.probe"]), 
#													 tumors[selected.probes.paired, , drop=FALSE], 
#													 normals[selected.probes.paired, , drop=FALSE], 
#													 regulation=regulation, 
#													 stddev, TRUE)
#		affected.samples.paired = list(summary=rbind(nonbody$summary, body$summary), samples=c(nonbody$samples, body$samples))
#  	}
  
#  gene.scores = list(unpaired=rep(0.0, length(genes)), paired=rep(0.0, length(genes))) 
#  names(gene.scores$unpaired) = genes
#  names(gene.scores$paired) = genes
#  samples = list(unpaired=list(), paired=list())
#  summary = list(unpaired=matrix(0, ncol=2, nrow=length(genes)), paired=matrix(0, ncol=2, nrow=length(genes)))
#  rownames(summary$unpaired) = genes
#  colnames(summary$unpaired) = c("absolute", "relative")
#  rownames(summary$paired) = genes
#  colnames(summary$paired) = c("absolute", "relative")
#  nTumors = ncol(tumors)
#  for (gene in genes) {
#    # Get all probes for that gene
#    geneProbes = gene2probe[[gene]]
#	allProbes = rownames(probe.annotation[geneProbes, , drop=FALSE])
#    
#    temp = length(geneProbes)
#    if (temp > 0) {
#	  gene.scores$unpaired[gene] = length(intersect(allProbes, selected.probes.unpaired)) / temp
#      gene.scores$paired[gene] = length(intersect(allProbes, selected.probes.paired)) / temp
#	  samples$unpaired[[gene]] = unique(unlist(affected.samples.unpaired$samples[intersect(allProbes, selected.probes.unpaired)]))
#	  samples$paired[[gene]] = unique(unlist(affected.samples.paired$samples[intersect(allProbes, selected.probes.paired)]))
#	  summary$unpaired[gene, "absolute"] = length(samples$unpaired[[gene]])
#	  summary$unpaired[gene, "relative"] = summary$unpaired[gene, "absolute"] / nTumors
#	  summary$paired[gene, "absolute"] = length(samples$paired[[gene]])
#	  summary$paired[gene, "relative"] = summary$paired[gene, "absolute"] / nTumors
#    }
#  }
#  
#  return(list(scores=gene.scores, diffs=mean.diff, cors=cors, wilcox=wilcox, corrected=wilcox.bh, samples=samples, summary=summary))
}


##' Calculates the ratio of significant probes per gene and the list of samples that are affected by changes in these probes.
##' @param tumors methylation matrix for tumor samples with probes in rows and samples in columns
##' @param normals methylation matrix for normal samples with probes in rows and samples in columns
##' @param meth.analysis list returned by doMethylationAnalysis()
##' @param gene2probe nested list mapping probes to genes and gene regions
##' @param probe2gene nested list mapping genes to probes and gene regions
##' @param genes character vector of genes to consider
##' @param wilcox.FDR significance cut-off for Wilcoxon FDR; default=0.05
##' @param cor.FDR significance cut-off for correlation FDR; default=0.05
##' @param diff.cutoff methylation difference needs to be smaller or greater than this cut-off to be considered significant; default=0.1
##' @param regulation either "down" or "up" for finding genes that are regulated in the corresponding direction; default="down"
##' @param gene.region boolean indicating whether gene region annotation should be taken into account; default=TRUE
##' @param stddev how many standard deviations does a sample have to be away from the mean to be considered affected; default=1
significantProbes = function(tumors,
				 			 normals,
							 meth.analysis,
							 gene2probe,
							 probe2gene,
							 genes, 
							 wilcox.FDR=0.05, 
							 cor.FDR=0.05,
							 diff.cutoff=0.1, 
							 regulation=c("down", "up"),
							 gene.region=TRUE,
							 stddev=1) { 
	regulation = match.arg(regulation)
	
	# Get the correct comparison function
	# If we want to find genes with higher methylation in tumors get the greaterThan function
	# If we want to find genes with lower methylation in tumors, get the smallerThan function
	compare = switch(regulation, down=greaterThan, up=smallerThan)
	
	# Correlation significance filter
	# Multiple testing correction
	meth.analysis$cors = cbind(meth.analysis$cors, cor.FDR=p.adjust(meth.analysis$cors$cor.p, method="BH"))
	significant.cors = meth.analysis$cors[which(meth.analysis$cors$cor.FDR <= cor.FDR), ]
	
	if (!gene.region) {
		significant.cors = cbind(significant.cors[which(significant.cors[, "cor"] < 0), ], region="rest")
	} else {
		tsc.body = data.frame()
		tsc.rest = data.frame()
		for (probe in unique(significant.cors[, "probe"])) {
			bodyGenes = probe2gene[[probe]][["Body"]]
			if (!is.null(bodyGenes)) {
				tsc.body = rbind(tsc, significant.cors[which(significant.cors[, "probe"] == probe & 
															 significant.cors[, "cor"] > 0 & 
															 significant.cors[, "gene"] %in% bodyGenes), ])
			}
			otherGenes = unlist(probe2gene[[probe]][setdiff(names(probe2gene[[probe]]), c("Body"))])
			if (!is.null(otherGenes)) {
				tsc.rest = rbind(tsc, significant.cors[which(significant.cors[, "probe"] == probe & 
													 		 significant.cors[, "cor"] < 0 & 
															 significant.cors[, "gene"] %in% otherGenes), ])
			}
		}
		significant.cors = rbind(cbind(tsc.body, region="body"), cbind(tsc.rest, region="rest"))
	}
	significant.probes = unique(significant.cors[, "probe"])
	
	# Calculate multiple testing correction for remaining probes
	meth.analysis$wilcox = cbind(wilcox.p=meth.analysis$wilcox, 
								 wilcox.FDR=p.adjust(meth.analysis$wilcox[significant.probes], method="BH")[names(meth.analysis$wilcox)])
	# Mean difference filter 
	significant.probes = intersect(significant.probes, names(meth.analysis$diffs)[compare(meth.analysis$diffs, diff.cutoff)])
	significant.probes = significant.probes[!is.na(significant.probes)]
	
	# Wilcoxon significance filter
	significant.probes = intersect(significant.probes, rownames(meth.analysis$wilcox)[meth.analysis$wilcox$wilcox.FDR <= wilcox.FDR])
	
	## Find affected samples per probe
	bodyprobes = intersect(significant.probes, significant.cors[which(significant.cors[, "region"] == "body"), "probe"])
	nonbodyprobes = intersect(significant.probes, significant.cors[which(significant.cors[, "region"] != "body"), "probe"])
	nonbody = countAffectedSamples(nonbodyprobes, 
								   tumors[nonbodyprobes, , drop=FALSE], 
								   normals[nonbodyprobes, , drop=FALSE], 
								   regulation=switch(regulation, up="down", down="up"), 
								   stddev, 
								   paired=FALSE)
	body = countAffectedSamples(bodyprobes, 
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
		allProbes = unlist(gene2probe[[gene]])
		# Get the ratio of (#significant probes with higher methylation) / (#all probes)
		gene.scores[gene] = length(intersect(significant.probes, geneProbes)) / length(allProbes)
		
		samples[[gene]] = union(unlist(body$samples[intersect(allProbes, names(body$samples))]),
								unlist(nonbody$samples[intersect(allProbes, names(nonbody$samples))]))
		summary[gene, "absolute"] = length(sampels[[gene]])
		summary[gene, "relative"] = summary[gene, "absolute"] / nTumors
	}
	
	list(scores=gene.scores, summary=summary, samples=samples)
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
summarizeMethylation = function(methylation, genes, threshold=0, score=c("atleast", "absolute")) {
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
