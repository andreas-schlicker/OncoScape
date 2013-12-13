# This file contains functions related to the scoring of prioritization results.
# 
# Author: Andreas Schlicker
###############################################################################


##' Comine scores from the given list into a matrix.
##' Missing scores are set to NA.
##' @param scores named list with all scores
##' @param genes vector with gene IDs; default: NULL (take all genes with any score)
##' @return matrix with scores; genes in rows and score categories in columns
##' @author Andreas Schlicker
combineScores = function(scores, genes=NULL) {
	if (is.null(genes)) {
		genes = c()
		for (n in names(scores)) {
			genes = union(genes, names(scores[[n]]))
		}
	}
	
	res = matrix(NA, ncol=length(scores), nrow=length(genes))
	colnames(res) = names(scores)
	rownames(res) = genes
	for (n in names(scores)) {
		res[intersect(genes, names(scores[[n]])), n] = scores[[n]][intersect(genes, names(scores[[n]]))]
	}
	
	res
}

##' Compute summarized scores for a given gene set. 
##' NA values are ignored. The summarize function has to accept a logical parameter "na.rm".
##' @param scoreMat matrix with score types in columns and genes in rows
##' @param geneset vector with gene IDs
##' @param summarize function for summarizing scores; default: mean
##' @return vector with summarized scores; vector of NA if scoreMat doesn't contain any gene
##' @author Andreas Schlicker
scoreSet = function(scoreMat, geneset, summarize=mean) {
	common = intersect(rownames(scoreMat), geneset)
	if (length(common) == 0) {
		res = rep(NA, times=ncol(scoreMat))
		names(res) = colnames(scoreMat)
	} else {
		res = apply(scoreMat[common, , drop=FALSE], 2, summarize, na.rm=TRUE)
	}
	
	res
}

##' Compute summarized scores for a list of gene sets.
##' NA values are ignored. The summarize function has to accept a logical parameter "na.rm".
##' @param scoreMat matrix with score types in columns and genes in rows
##' @param genesets named list with vectors containing gene IDs
##' @param summarize function for summarizing scores; default: mean
##' @return matrix with summarized scores, gene sets in rows and score types in rows
##' @author Andreas Schlicker
scoreGenesets = function(scoreMat, genesets, summarize=mean) {
	if (is.null(names(genesets))) {
		names(genesets) = as.character(1:length(genesets))
	}
	
	sumScores = lapply(genesets, function(x) { scoreSet(scoreMat, x, summarize)})
	names(sumScores) = names(genesets)
	
	scoreNames = colnames(scoreMat)
	res = matrix(NA, ncol=ncol(scoreMat), nrow=length(genesets))
	colnames(res) = scoreNames
	rownames(res) = names(genesets)
	for (n in names(sumScores)) {
		res[n, ] = sumScores[[n]][scoreNames]
	}
	
	res
}
