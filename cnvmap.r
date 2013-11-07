##' Queries Biomart to obtain the chromosomal location for each input gene.
##' @param genes character vector with HGNC symbols
##' @param warningsfile file to save warning messages in; default: "warnings.txt".
##' If the file doesn't exist, it will be created. Otherwise, new warnings are appended.
##' @param genome genome version to be used; default: hg18
##' @return Returns a data.frame with the chromosome, start and end positions
##' @author Andreas Schlicker
getGeneLocs = function(genes, warningsfile="warnings.txt", genome="hg18") {
	require(biomaRt) || stop("Can't load package \"biomaRt\"!")
	
	# List all Biomarts available through the corresponding archive
	# listMarts(host="may2009.archive.ensembl.org",path="/biomart/martservice")
	
	# Use the latest Biomart that is based on NCBI36/hg18
	mart = useMart(host="may2009.archive.ensembl.org", path="/biomart/martservice", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
	if (genome == "hg19") {
		mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
	}
	
	# Find out the gene locations
	res = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters=c("hgnc_symbol"), values=genes, mart=mart)
	res = subset(res, chromosome_name %in% c(1:22, "X", "Y"))
	
	wf = NULL
	
	if (length(unique(res$hgnc_symbol)) < length(res$hgnc_symbol)) {
		if (is.null(wf)) {
			wf = file(warningsfile, "a")
		}
		cat(paste("Found several genomic locations for genes: ", res$hgnc_symbol[which(duplicated(res$hgnc_symbol))], "\n", sep=""), file=wf)
		cat("Using the first location for each gene!\n", file=wf)
		res = res[!duplicated(res$hgnc_symbol), ]
	}
	
	if (!is.null(wf)) {
		close(wf)
	}
	
	rownames(res) = res[, "hgnc_symbol"]
	if (length(which(res[, "chromosome_name"] == "X")) > 0) {
		res[which(res[, "chromosome_name"] == "X"), "chromosome_name"] = 23
	}
	if (length(which(res[, "chromosome_name"] == "Y")) > 0) {
		res[which(res[, "chromosome_name"] == "Y"), "chromosome_name"] = 24
	}
	
	res
}

##' Find the correct copy number segment and the associated copy number.
##' @param start start location
##' @param stop end location
##' @param segments matrix of start (column 1) and end (column 2) positions, and the
##' copy number value (column 3) of the segments; all segments have to be on the
##' same chromosome as the position of interest
##' If a gene spans several segments, the mean copy number of all these segments
##' is returned.
##' @return named list with two elements: cnvalue: copy number value; case: integer 
##' code; 1=gene is not part of any segment; 2=gene is part of exactly one segment;
##' 3=start of gene lies between segments; 4=end of gene lies between segments;
##' 5=gene spans several segments
##' @author Andreas Schlicker
getCopyNumberValue = function(start, stop, segments) {
	if (!is.matrix(segments)) {
		# Segments is a vector with only one segment, so add one artificial to make 
		# it a matrix. 
		segments = rbind(c(segments[1], 0, 0, 0), segments)
	}
	
	# Find the start segment
	# Look for the largest segment start that is still smaller than the gene start
	startSeg = 0
	startSeg = max(which(segments[, 2] <= start))
	# If the end of that segment is larger than the start, the gene doesn't touch
	# that segment.
	startWithin = TRUE
	if (is.infinite(startSeg) || segments[startSeg, 3] < start)
		startWithin = FALSE
	
	# and the end segment
	# Look for the smallest segment end that is still bigger than the gene end
	endSeg = 0
	endSeg = min(which(segments[, 3] >= stop))
	# If the start of that segment is larger than the end, the gene doesn't touch
	# that segment.
	endWithin = TRUE
	if (is.infinite(endSeg) || segments[endSeg, 2] > stop)
		endWithin = FALSE
	
	
	if (is.infinite(startSeg) && !is.infinite(endSeg) || (startSeg > endSeg)) {
		# If no start segment but an end segment was found, the gene begins before the first one.
		# Use the first segment in this case
		startSeg = 1
	} else if (!is.infinite(startSeg) && is.infinite(endSeg)) {
		# If no end segment but a start segment was found, the gene ends behind the last one.
		# Use the last segment in this case
		endSeg = nrow(segments)
	}
	
	case = 0
	if (is.infinite(startSeg) && is.infinite(endSeg)) {
		# The gene doesn't touch any segments
		cnvalue = NA
		case = 1
	} else if (startSeg == endSeg) {
		# The gene lies withing one segment
		cnvalue = segments[startSeg, 4]
		case = 2
	} else if (!startWithin) {
		# The gene's start lies between segments
		# -> ignore the start segment and average the remaining
		cnvalue = mean(segments[(startSeg+1):endSeg, 4])
		case = 3
	} else if (!endWithin) {
		# The gene's end lies between segments
		# -> ignore the end segment and average the remaining
		cnvalue = mean(segments[startSeg:(endSeg-1), 4])
		case = 4
	} else {
		# The gene spans different segments
		# Take the average of all copy number values of these segments
		cnvalue = mean(segments[startSeg:endSeg, 4])
		case = 5
	}
	
	list(cnvalue=as.double(cnvalue), case=case)
}

##' Maps all input genes to the copy number value of the corresponding segment.
##' @param genes character vector containing HGNC symbols
##' @param cnData matrix with segmented copy number variation data; formatted
##' as CBS output
##' @param cases boolean; if set to TRUE, an additional matrix containing 
##' information on how each CN value was determined is returned; default: FALSE
##' @param genome version of the genome; default: "hg18"
##' @return named list with up to two elements; cn.values: gene-by-sample copy
##' number matrix; cn.cases: gene-by-sample conversion case matrix
##' @author Andreas Schlicker
mapGenes2CN = function (genes, cnData, cases=FALSE, genome="hg18") {
	# Map all genes to their chromosomal location using biomaRt
	genes2loc = getGeneLocs(genes, genome=genome)
	# List of samples
	allSamples = unique(cnData[, 1])
	
	# Matrix that will contain the copy number for each gene in each sample
	cn.matrix = matrix(NA, nrow=length(genes), ncol=length(allSamples))
	rownames(cn.matrix) = genes
	colnames(cn.matrix) = allSamples
	
	# If the case matrix is to be reported
	if (cases) {
		cases = matrix(-1, nrow=length(genes), ncol=length(allSamples))
		rownames(cases) = genes
		colnames(cases) = allSamples
	}
	
	# Go through all samples
	for (samp in allSamples) {
		# Get all segments with their copy number for the current sample
		samp.segs = as.matrix(cnData[which(cnData[, 1] == samp), c(2, 3, 4, 6)])
		# Go through the chromosomes
		for (chrom in unique(genes2loc[, "chromosome_name"])) {
			# Get the segments for that particular chromosome
			segs = samp.segs[which(samp.segs[, 1] == chrom), ]
			# Go through all genes on that chromosome
			for (gene in rownames(genes2loc[which(genes2loc[, "chromosome_name"] == chrom), ])) {
				# Get the copy numberk
				res = getCopyNumberValue(as.integer(genes2loc[gene, "start_position"]), as.integer(genes2loc[gene, "end_position"]), segs)
				cn.matrix[gene, samp] = res[["cnvalue"]]
				# Remember the case if that information is wanted
				if (cases) {
					cases[gene, samp] = res[["case"]]
				}
			}
		}
	}
	
	res = list(cn.values=cn.matrix)
	if (cases) {
		res[["cn.cases"]] = cases
	}
	
	res
}

##' Calculates correlation between copy number and expression of the corresponding genes.
##' Correlation is calculated using all samples contained in both data matrices
##' @param cnv.data matrix with genes in rows and samples in columns
##' @param exprs.data matrix with genes in rows and samples in columns
##' @param genes character vector with genes to test
##' @return data.frame with correlation values and p-values
corCnvExprs = function(cnv.data, exprs.data, genes) {
	# All genes contained in the expression data
	common.genes = intersect(genes, intersect(rownames(cnv.data), rownames(exprs.data)))
	commonSamples = intersect(colnames(cnv.data), colnames(exprs.data))
	
	cnv.data = cnv.data[, commonSamples, drop=FALSE]
	exprs.data = exprs.data[, commonSamples, drop=FALSE]
	
	cors = data.frame(gene=character(length(genes)), cor=rep(NA, times=length(genes)), 
				      cor.p=rep(NA, times=length(genes)), stringsAsFactors=FALSE)
	rownames(cors) = genes
	for (x in common.genes) {
		tempCor = cor.test(exprs.data[x, ], 
						   cnv.data[x, ], 
						   method="spearman", 
						   use="pairwise.complete.obs",
						   exact=FALSE)
		cors[x, ] = c(x, tempCor$estimate, tempCor$p.value)
	}
	
	cors[, 2] = as.numeric(cors[, 2])
	cors[, 3] = as.numeric(cors[, 3])
	
	cors
}

##' Run copy number alteration analysis for each gene. 
##' @param tumors matrix with tumor copy number matrix, genes in rows and samples in columns
##' @param normals matrix with normal sample copy number matrix, genes in rows and samples in columns
##' @param exprs expression data matrix, genes in rows and samples in columns
##' @param genes character vector of genes; default: NULL (test all genes in both tumors and normals)
##' @param samples vector with sample names to use for the analysis. If this is NULL, all samples will
##' be used; default: NULL
##' @param paired.wilcox boolean indicating whether doing paired or unpaired analysis; default: TRUE
##' @return named list with gene scores, correlation results, Wilcoxon test results and difference in mean
##' @author Andreas Schlicker
doCnvAnalysis = function(tumors, 
						 normals, 
						 exprs, 
						 genes=NULL,
						 samples=NULL,
						 paired.wilcox=TRUE) {
  
  	if (is.null(samples)) {
    	tumor.samples = colnames(tumors)
    	normal.samples = colnames(normals)
  	} else {
	    tumor.samples = intersect(samples, colnames(tumors))
    	normal.samples = intersect(samples, colnames(normals))
  	}

	selected.genes = intersect(rownames(tumors), rownames(normals))
	if (!is.null(genes)) {
		selected.genes = intersect(genes, selected.genes)
	}	
  	if (length(selected.genes) == 0) {
	  	warning("No gene of interest is contained in copy number data of both tumors and normals!")
	  	return(list())
  	}
  
	tumors = tumors[selected.genes, tumor.samples, drop=FALSE]
  	normals = normals[selected.genes, normal.samples, drop=FALSE]
  
  	# Calculate the difference in copy number between tumors and normals
  	# tumor - normal
  	# --> diff < 0 implies that the gene has higher mean copy number values in normals than in tumors
  	# --> diff > 0 implies that the gene has lower mean copy number values in normals than in tumors
  	mean.diff = meanDiff(tumors, normals)
  	# Calculate correlation between copy number and expression
  	cors = corCnvExprs(tumors, exprs, selected.genes)
  	# Run Wilcoxon tests
  	wilcox = doWilcox(tumors, normals, paired.wilcox)
  
	list(diffs=mean.diff, cors=cors, wilcox=wilcox)
}

##' Score genes and find out which samples are affected by corresponding aberrations.
##' @param tumors matrix with tumor copy number matrix, genes in rows and samples in columns
##' @param normals matrix with normal sample copy number matrix, genes in rows and samples in columns
##' @param cnv.analysis list returned by doCnvAnalysis()
##' @param genes character vector of genes; default: NULL (test all genes in both tumors and normals)
##' @param wilcox.FDR significance cut-off for Wilcoxon FDR; default=0.05
##' @param cor.FDR significance cut-off for correlation FDR; default=0.05
##' @param diff.cutoff copy number difference needs to be smaller or greater than this cut-off to be considered significant; default=-0.1
##' @param regulation either "down" or "up" for finding genes that are regulated in the corresponding direction; default="down"
##' @param stddev how many standard deviations does a sample have to be away from the mean to be considered affected; default=1
##' @param paired boolean indicating whether paired or unpaired analysis was performed; default: TRUE
##' @return named list with scores for genes ("scores"), number of affected samples ("summary") and the lists of affected samples ("samples")
##' @author Andreas Schlicker
summarizeCnv = function(tumors,
						normals,
						cnv.analysis, 
						genes=NULL, 
						wilcox.FDR=0.05, 
						cor.FDR=0.05, 
						diff.cutoff=-0.1, 
						regulation=c("down", "up"),
						stddev=1,
						paired=TRUE) {
	regulation = match.arg(regulation)
					
	# Get the correct comparison function
	# If we want to find genes with higher copy number in tumors get the greaterThan function
	# If we want to find genes with lower copy number in tumors, get the smallerThan function
	compare = switch(regulation, down=smallerThan, up=greaterThan)
	
	significant.genes = intersect(rownames(tumors), rownames(normals))
	if (!is.null(genes)) {
		significant.genes = intersect(genes, significant.genes)
	}
	cnv.analysis$cors = cnv.analysis$cors[which(rownames(cnv.analysis$cors) %in% significant.genes), ]
	cnv.analysis$wilcox = cnv.analysis$wilcox[which(names(cnv.analysis$wilcox) %in% significant.genes)]
	cnv.analysis$diffs = cnv.analysis$diffs[which(names(cnv.analysis$diffs) %in% significant.genes)]
	
	## Integrate results
	# Correlation significance filter
	cnv.analysis$cors = cbind(cnv.analysis$cors, cor.FDR=p.adjust(cnv.analysis$cors[, "cor.p"], method="BH"))
	significant.genes = rownames(cnv.analysis$cors)[which(cnv.analysis$cors[, "cor.FDR"] <= cor.FDR & cnv.analysis$cors[, "cor"] > 0)]
	
	# Difference filter
	cnv.analysis$wilcox = cbind(wilcox.p=cnv.analysis$wilcox, 
								wilcox.FDR=p.adjust(cnv.analysis$wilcox[significant.genes], method="BH")[names(cnv.analysis$wilcox)])
	significant.genes = names(which(compare(cnv.analysis$diffs[significant.genes], diff.cutoff)))
	significant.genes = names(which(cnv.analysis$wilcox[significant.genes, "wilcox.FDR"] <= wilcox.FDR))
	
	gene.scores = rep(0, length(genes))
	names(gene.scores) = genes
	gene.scores[significant.genes] = 1
	
	affected.sampels = countAffectedSamples(significant.genes, tumors, normals, regulation, 1, TRUE)
	
	list(scores=gene.scores, summary=affected.samples$summary, samples=affected.samples$samples)
}
