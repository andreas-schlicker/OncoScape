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
			function(cancType) { list(prioritize.combined=results1[[cancType]]$prioritize.combined - 
														  results2[[cancType]]$prioritize.combined[rownames(results1[[cancType]]$prioritize.combined), ]) })
	names(res) = cancTypes
	
	res
}

##' Returns a vector with scores summed across all cancer types.
##' The vector is sorted in ascending order.
##' @param scoreDf score data frame
##' @param score name of the score to use
##' @return named vector with sum of changes and gene names
##' @author Andreas Schlicker
alteredGenes = function(scoreDf, score="Combined") {
	scoreDf = scoreDf[which(scoreDf[, "score.type"] == score & scoreDf[, "score"] != 0), ]
	sort(unlist(lapply(split(scoreDf$score, scoreDf$gene), sum, na.rm=TRUE)))
}

##' Generates four plots summarizing the changes in a particular score.
##' @param results1 list with baseline results
##' @param results2 list with changed results
##' @param score named list with the score to be summarized;
##' The list element has to match a column name in the results matrices.
##' The name of the element will be used as label.
##' Combined.affected is a special score, which is translated to og.affected.rel - ts.affected.rel
##' @return named list with four plots: (un-)faceted confusion heatmaps and (un-)faceted difference barplots
##' @author Andreas Schlicker
compareScore = function(results1, results2, score) {
	res1.hdt = heatmapDataframe(results1, score[1])
	res2.hdt = heatmapDataframe(results2, score[1])
	combined = data.frame()
	for (cancType in as.character(unique(res1.hdt[, "cancer"]))) {
		temp1 = res1.hdt[which(res1.hdt[, "cancer"] == cancType), ]
		temp2 = res2.hdt[which(res2.hdt[, "cancer"] == cancType), ]
		combined = rbind(combined, 
						 cbind(temp1,
							   score2=temp2[match(temp1[, "gene"], temp2[, "gene"]), "score"], 
							   score.diff=temp2[match(temp1[, "gene"], temp2[, "gene"]), "score"]))
	}
	
	p1 = confusionHeatmap(conHeatDf(combined), xlab="baseline", ylab="altered score", title="")
	p2 = barplot(combined)
	p3 = confusionHeatmap(conHeatDf(combined), facet=NULL, xlab="baseline", ylab="altered score", title="")
	p4 = barplot(combined, facet=NULL)
	
	list(confusion=p1, difference=p2, confusionSum=p3, differenceSum=p4)
}
