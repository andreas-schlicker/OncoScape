# Functions for performing differential expression analysis using edgeR on
# sequencing data. 
#
# Author: Andreas Schlicker
###############################################################################

##' Creates a DGEList object used by other functions.
##' Rownames of the count.table matrix will be used as rownames for the DGEList object
##' @param count.table numeric read count matrix
##' @param feature.ann annotation for features; default: NULL
##' @param remove.zeros boolean, whether features with only zero counts should
##' be removed
##' @return the DGEList object
##' @author Andreas Schlicker
createDGEList = function(count.table, feature.ann=NULL, remove.zeros=FALSE) {
	if (!require(edgeR)) {
		stop("Could not load required package \"edgeR\"!")
	}
	dgel = DGEList(counts=count.table, genes=gene.ann, remove.zeros=remove.zeros)
	rownames(dgel) = rownames(count.table)
	
	dgel
}

##' Filters out features that are not expressed in given number of samples
##' @param dgel the DGEList object
##' @param count.cutoff count per million cutoff value that is applied to each sample
##' @param sample.cutoff at least sample.cutoff many samples need to meet the 
##' count.cutoff for a feature to be kept; default: 0.1 (= 10% of samples)
##' @param relative boolean, whether sample.cutoff is to be interpreted as a relative
##' fraction of samples or an absolute number; default: TRUE
##' @return filtered DGEList object
##' @author Andreas Schlicker
filterDGEList = function(dgel, count.cutoff, sample.cutoff=0.1, relative=TRUE) {
	# Convert absolute numer of samples to relative quantity
	if (!relative) {
		sample.cutoff = sample.cutoff / ncol(dgel)
	}
		
	keep = (rowSums(cpm(dgel) > count.cutoff) / ncol(dgel)) >= sample.cutoff
	dgel = dgel[keep, ]
	# Correct library size after filtering
	dgel$samples$lib.size = colSums(dgel$counts)
	
	dgel
}

##' Creates a design matrix.
##' @param sample.names character vector with the sample names; will be used as
##' rownames for the design matrix. Has to contain all tumor samples and then all
##' normal samples. Each entry in this vector has to be unique.
##' @param tumors character vector with the names of the tumor samples 
##' @param normals character vector with the names of the normal samples
##' @only.paired boolean, indicates whether only samples with tumor and normal
##' should be included (default: TRUE). If TRUE, normal and tumor of the 
##' same sample need to have the same name
##' @return the design matrix
##' @author Andreas Schlicker
designMatrix = function(sample.names, tumors, normals, only.paired=TRUE) {
	patients = tumors
	if (only.paired) {
		patients = intersect(tumors, normals)
	}
	tissue = factor(c(rep("T", times=length(patients)), 
					  rep("N", times=length(normals))))
	# Add the normal samples
	patients = factor(c(patients, normals))
	design = model.matrix(~patients+tissue)
	rownames(design) = sample.names
	
	design
}

##' Runs the different dispertion estimation steps.
##' @param dgel DGEList object
##' @param design the design matrix
##' @return the updated DGEList object
##' @author Andreas Schlicker
estimateDispertion = function(dgel, design) { 
	# Estimate the common dispersion
	dgel =  estimateGLMCommonDisp(dgel, design, verbose=FALSE)
	# gene-wise dispersion values
	dgel <- estimateGLMTrendedDisp(dgel, design)
	dgel <- estimateGLMTagwiseDisp(dgel, design)
	dgel
}

##' Runs the differential expression analysis.
##' @param dgel DGEList object
##' @param design design matrix
##' @return data frame with the results of the analysis
##' @author Andreas Schlicker
diffExpr = function(dgel, design) {
	fit = glmFit(dgel, design)
	
	topTags(lrt, n=nrow(lrt))[[1]]
}

##' Compares expression data of tumor and normal samples.
##' The function will run (un-)paired Wilcoxon tests. Pairing of samples
##' is determined by matching sample names in the tumors and normals matrices.
##' The function uses doWilcox() as defined for the analysis of methylation data!
##' @param tumors expression matrix with samples in the columns and genes in the rows
##' @param normals expression matrix with samples in the columns and genes in the rows
##' @param paired boolean, whether paired or unpaired test is to be performed
##' @return names character vector with results from tests
##' @author Andreas Schlicker
runExprComp = function(tumors, normals, paired=TRUE) {
	if (!is.matrix(tumors)) {
		tumors = matrixFromVector(tumors)
	}
	if (!is.matrix(normals)) {
		normals = matrixFromVector(normals)
	}
	
	if (paired) {
		# Samples from tumors that have a matched normal
		matchedsamples = intersect(colnames(tumors), colnames(normals))
		if (length(matchedsamples) > 0) {
			tumors = tumors[, matchedsamples, drop=FALSE]
			normals = normals[, matchedsamples, drop=FALSE]
		} else {
			paired = FALSE
			warning("No paired expression samples found. Performing unpaired analysis!")
		}
	}
	
	# Rename normal samples to reflect the state
	colnames(normals) = paste(colnames(normals), "normal", sep="_")
	
	sample.groups = c(rep(1, times=ncol(tumors)), rep(2, times=ncol(normals)))
	names(sample.groups) = c(colnames(tumors), colnames(normals))
	
	# Run Wilcoxon tests
	if (paired) {
		wilcox.p = doWilcox(cbind(tumors, normals), matchedSamples=matchedsamples)
	} else {
		wilcox.p = doWilcox(cbind(tumors, normals), groups=sample.groups)
	}
	
	data.frame(tumor.exprs=apply(tumors, 1, mean, na.rm=TRUE),
			   normal.exprs=apply(normals, 1, mean, na.rm=TRUE),
			   p.value=wilcox.p,
			   BH=p.adjust(wilcox.p, method="BH"))
}

##' Tests for differential expression between tumors and normals.
##' @param genes character vector of gene symbols
##' @param exprs.res data.frame with expression differences as returned by runExprComp()
##' @param regulation either "down" or "up", whether down- or upregulation in tumors
##' should be scored
##' @param wilcox.cutoff FDR threshold below which a difference is considered to be
##' significant
##' @return names list with the results
##' @author Andreas Schlicker
doExprAnalysis = function(genes, exprs.res, regulation=c("down", "up"), wilcox.cutoff=0.05) {
	regulation = match.arg(regulation)
	
	# Get the correct comparison function
	# If we want to find genes with greater expression in tumors get the greaterThan function
	# If we want to find genes with lower expression in tumors, get the smallerThan function
	compare = switch(regulation, down=smallerThan, up=greaterThan)
	
	gene.scores = rep(0, length(genes))
	names(gene.scores) = genes
	hits = apply(exprs.res[intersect(rownames(exprs.res), genes), ], 
				 1, function(x) { x["BH"] < wilcox.cutoff & compare(x["tumor.exprs"], x["normal.exprs"]) })
	gene.scores[names(hits)[hits]] = 1
	
	gene.scores
}
