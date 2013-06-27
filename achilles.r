##' Fetch the values for the given gene and cell lines
##' from the matrix
##' @param mat the matrix
##' @param gene the gene symbol
##' @param cls the cell lines to look for
##' @return a vector with all phenotype values; NA if the 
##' gene is not contained in the data set
##' @export 
##' @author Andreas Schlicker
fetchValues = function(mat, gene, cls) {
	res = NA
	if (gene %in% rownames(mat)) {
		res = mat[gene, cls]
	}
	
	res
}

##' Fetches the phenotype values of the given gene for
##' the given cell lines.
##' @param gene the gene symbol
##' @param cls names of the cell lines
##' @return a vector with all phenotype values; NA if the 
##' gene is not contained in the data set
##' @export 
##' @author Andreas Schlicker
getPhenoValues = function(gene, cls) {
	fetchValues(achilles, gene, cls) 
}

##' Fetches the rank of the given gene for 
##' the given cell lines.
##' @param gene the gene symbol
##' @param cls the cell line names
##' @param relative boolean, if true, the relative rank will
##' be returned, else the absolute rank (default)
##' @return a vector with all ranks; NA if the 
##' gene is not contained in the data set
##' @export 
##' @author Andreas Schlicker
getGeneRanks = function(gene, cls, relative=FALSE) {
	if (relative) {
		res = fetchValues(achilles.rank.rel, gene, cls)
	} else {
		res = fetchValues(achilles.rank, gene, cls)
	}
	
	res
}

##' Fetches the absolute or relative average rank for the given gene
##' across all cell lines.
##' @param gene the gene symbol
##' @param relative boolean, if true, the relative rank will
##' be returned, else the absolute rank (default)
##' @return the absolute or relative rank for the gene or NA if the 
##' gene was not found
##' @export 
##' @author Andreas Schlicker
getGeneAvgRank = function(gene, relative=FALSE) {
	column = "avg.ranK"
	if (relative) {
		column = "avg.rank.rel"
	}
	
	res = NA
	if (gene %in% rownames(achilles.rank.avg)) {
		res = achilles.rank.avg[gene, column]
	}
	
	res
}

##' Get the names of cell lines for the given tissue.
##' @param tissue the tissue
##' @return a vector with the corresponding cell line
##' names; all cell line names if no cell line was found
##' @export 
##' @author Andreas Schlicker
getCellLines = function(tissue) {
	cls = rownames(subset(achilles.clann, Site_primary==tissue))
	if (length(cls) == 0) {
		cls = rownames(achilles.clann)
	}
	
	cls
}

##' Wanna know which tissue types are available?
##' Use this function.
##' @return a character vector with unique tissues
##' @export
##' @author Andreas Schlicker
getTissues = function() {
	unique(achilles.clann[, "Site_primary"])
}

##' Loads the data for Project Achilles
##' @param the directory
##' @export 
##' @author Andreas Schlicker
loadAchillesData = function(directory) {
	load(paste(directory, "achilles_clann.rdata", sep="/"), envir=.GlobalEnv)
	load(paste(directory, "achilles_rank_avg.rdata", sep="/"), envir=.GlobalEnv)
	load(paste(directory, "achilles_rank.rdata", sep="/"), envir=.GlobalEnv)
	load(paste(directory, "achilles_rank_rel.rdata", sep="/"), envir=.GlobalEnv)
	load(paste(directory, "achilles.rdata", sep="/"), envir=.GlobalEnv)
}

##' Runs the analysis using data from Project Achilles.
##' @param genes a character vector with gene symbols to look at
##' @param score use the "rank", "avgrank" or the "phenotype" 
##' scores of the genes for the analysis. "rank" refers to the
##' absolute rank of the gene. "avgrank" refers to the mean rank the
##' gene achieved across all cell lines in the analysis. "phenotype"
##' refers to the phenotype score provided by Project Achilles. Smaller
##' ranks and scores indicate higher survival dependence on the gene. 
##' default: "rank"
##' @param tissue the tissue type of cell lines to include.
##' Use getTissues() to get a list of available tissue types.
##' If there are no cell lines for the given tissue, all cell lines
##' will be included in the analysis. default: "all"
##' @param relative boolean; indicates whether relative ranks should
##' be used for the analysis. Is ignored if score=="phenotype". default: TRUE
##' @export 
##' @author Andreas Schlicker
doAchillesAnalysis = function(genes, score=c("rank", "avgrank", "phenotype"), tissue="all", relative=TRUE) {
	# Make sure the score argument is valid
	score = match.arg(score)
	# Get the cell lines for the tissue
	cls = getCellLines(tissue)
	
	if (score == "rank") {
		scores = lapply(genes, getGeneRanks, cls=cls, relative=relative)
	} else if (score == "avgrank") {
		cls = getCellLines(tissue)
		if (length(cls) == ncol(achilles)) {
			scores = lapply(genes, getGeneAvgRanks, relative=relative)
		} else {
			scores = lapply(genes, getGeneRanks, cls=cls, relative=relative)
			scores = lapply(scores, mean)
		}
	} else if (score == "phenotype") {
		scores = lapply(genes, getPhenoValues, cls=cls)
	}
	
	unlist(scores)
}

##' Summarizes the data from Project Achilles. Basically compares
##' whether values for single genes are less/greater than the given 
##' threshold. 
##' @param achil.scores a vector with phenotype scores or ranks; typically
##' obtained from doAchillesAnalysis()
##' @param threshold the threshold to compare to
##' @param alternative "less" than threshold or "greater" than treshold; default: "less"
##' @return a named vector with 1 if the test for the gene was true, 0 otherwise
##' @export
##' @author Andreas Schlicker
summarizeAchilles = function(achil.scores, threshold, alternative=c("less", "greater")) {
	alternative = match.arg(alternative)
  
  # Get the correct comparison function
  # If we want to find genes with positive effect on cell line growth, get the greaterThan function
  # If we want to find genes with positive effect on cell line growth, get the smallerThan function
  compare = switch(alternative, greater=greaterThan, less=smallerThan)
  
  sapply(achil.scores, function(x) { as.integer(compare(x, treshold)) })
}
