# Module that provides functions for comparing results from different prioritization
# runs.
# 
# Author: Andreas Schlicker
###############################################################################

##' Calculates the difference between the result matrices for each cancer type
##' contained in both result sets. The first one is used as base result.
##' @param results1 first result set
##' @param results2 second result set
##' @return named list of matrices with the differences
diffMat = function(results1, results2) {
	cancTypes = intersect(names(results1), names(results2))
	res = lapply(cancTypes, 
			function(cancType) { list(prioritize.combined=results[[cancType]]$prioritize.combined - results[[cancType]]$prioritize.combined) })
	names(res) = cancTypes
	
	res
}

##' Returns a vector with genes that have a non-zero score.
##' @param scoreDf score data frame
##' @param score name of the score to use
##' @return vector with gene names
##' @author Andreas Schlicker
alteredGenes = function(scoreDf, score="Combined") {
	as.character(unique(scoreDf[which(scoreDf[, "score.type"] == score & scoreDf[, "score"] != 0), "gene"]))
}