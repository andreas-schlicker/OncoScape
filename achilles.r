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
	column = "avg.rank"
	if (relative) {
		column = "avg.rel.rank"
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
	cls = rownames(subset(achilles.clann, Type==tissue))
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
	sort(unique(achilles.clann[, "Type"]))
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
##' The function either gets all the ranks or phenotypes for the 
##' selected cell lines. Then uses the function in the summarize 
##' argument to calculate one final value for each gene. 
##' @param genes a character vector with gene symbols to look at; default: NULL (all genes in matrix achilles)
##' @param score use the "rank" or the "phenotype" scores of the 
##' genes for the analysis. "rank" refers to the rank of the gene
##' for each cell line. "phenotype" refers to the phenotype score 
##' provided by Project Achilles. Smaller ranks and scores indicate 
##' higher survival dependence on the gene. default: "rank"
##' @param summarize function to summarize values for the different
##' cell lines. In order to count the number of values that are greater
##' or less than a certain cutoff, gtCutoff(), gtCutoffPercent(), 
##' ltCutoff() and ltCutoffPercent() can be used. default: median
##' @param cls character vector with cell line names that should be
##' included in the analysis. If null, use cell lines defined by the
##' "tissue" argument. 
##' @param tissue the tissue type of cell lines to include.
##' Use getTissues() to get a list of available tissue types.
##' If there are no cell lines for the given tissue, all cell lines
##' will be included in the analysis. default: "all"
##' @param relative boolean; indicates whether relative ranks should
##' be used for the analysis. Is ignored if score=="phenotype". default: TRUE
##' @return vector with the scores for each gene
##' @export 
##' @author Andreas Schlicker
doAchillesAnalysis = function(genes=NULL, score=c("rank", "phenotype"), summarize=median, cls=NULL, tissue="all", relative=TRUE) {
	# Make sure the score argument is valid
	score = match.arg(score)
	
	if (is.null(cls)) {
		# Get the cell lines for the tissue
		cls = getCellLines(tissue)
	}
	
	if (is.null(genes)) {
		genes = rownames(achilles)
	}
	
	if (score == "rank") {
		scores = lapply(genes, getGeneRanks, cls=cls, relative=relative)
	} else if (score == "phenotype") {
		scores = lapply(genes, getPhenoValues, cls=cls)
	}
	
	scores = unlist(lapply(scores, summarize))
	names(scores) = genes
	
	scores
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
  
  sapply(achil.scores, function(x) { as.integer(compare(x, threshold)) })
}
