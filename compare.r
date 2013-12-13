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
			function(cancType) { list(prioritize.combined=results1[[cancType]]$prioritize.combined - results2[[cancType]]$prioritize.combined) })
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

compareScore = function(results1, results2, score) {
	res1.hdt = heatmapDataframe(results1, score)
	res2.hdt = heatmapDataframe(results2, score)
	combined = data.frame()
	for (cancType in unique(res1.hdt[, "cancer"])) {
		temp1 = res1.hdt[which(res1.hdt[, "cancer"] == cancType), ]
		temp2 = res2.hdt[which(res2.hdt[, "cancer"] == cancType), ]
		combined = rbind(combined, 
						 cbind(temp1,
							   score2=temp2[match(temp1[, "gene"], temp2[, "gene"]), "score"], 
							   score.diff=temp2[match(temp1[, "gene"], temp2[, "gene"]), "score"]))
	}
	
	p1 = confusionHeatmap(conHeatDf(combined))
	p2 = barplot(combined)
	p3 = confusionHeatmap(conHeatDf(combined), facet=NULL)
	p4 = barplot(combined, facet=NULL)
	
	list(confusion=p1, difference=p2, confusionSum=p3, differenceSum=p4)
}