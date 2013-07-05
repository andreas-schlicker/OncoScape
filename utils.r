# Tests whether the given gene is expressed in at least in the given percentage of samples
# exprs.mat: expression matrix with genes in rows and samples in columns
# gene: gene to test
# exprs.threshold: threshold above which a gene is considered to be expressed; default: 0;
# only sensible for log-ratio expression data
# cutoff: percent samples that need to show expression for the gene to be considered expressed; 
# default: 0.5, i.e. considered expressed if expression in 50% of samples is > exprs.threshold
# returns: TRUE/FALSE
geneExpressed = function(exprs.mat, gene, exprs.threshold=0, cutoff=0.5) {
	(sum(sapply(exprs.mat[gene, ], function(x) { x > exprs.threshold })) / ncol(exprs.mat)) >= cutoff
}

# Convenience function to test a number of genes whether they are expressed or not
# exprs.mat: expression matrix with genes in rows and samples in columns
# gene: gene to test
# exprs.threshold: threshold above which a gene is considered to be expressed; default: 0;
# only sensible for log-ratio expression data
# cutoff: percent samples that need to show expression for the gene to be considered expressed; 
# default: 0.5, i.e. considered expressed if expression in 50% of samples is > exprs.threshold
# returns: a boolean vector
genesExpressed = function(exprs.mat, genes, exprs.threshold=0, cutoff=0.5) { 
	sapply(genes, function(gene) { geneExpressed(exprs.mat, gene, exprs.threshold, cutoff) })
}

# Perform Wilcoxon tests on the given matrix.
# If the matchedSamples argument is defined, a paired test is performed. In this
# case, the first length(matchedSamples) number of columns need to contain group
# 1 and the remaining columns samples in group 2. Samples in groups 1
# and 2 need to be in matched order.
# If a non-paired test is to be performed, groups should be a named vector of
# 1 and 2 indicating which samples belong to group 1 and 2, respectively.
doWilcox = function(inpMat, matchedSamples=NULL, groups=NULL) {
  # Paired Wilcoxon test?
  paired = !is.null(matchedSamples)
 
  # Get the two groups of samples
  if (paired) {
    group1 = 1:length(matchedSamples) 
    group2 = (length(matchedSamples)+1):ncol(inpMat)
  } else {
    group1 = which(colnames(inpMat) %in% names(groups[groups == 1]))
    group2 = which(colnames(inpMat) %in% names(groups[groups == 2]))
  }
  
  wilcox.p = apply(inpMat, 1, function(x) { tryCatch(wilcox.test(x[group1], x[group2], paired=paired, exact=FALSE)$p.value, error = function(e) NA) })
  
  return(wilcox.p)
}

# Perform Bartlett's test on each row of the given matrix.
# Groups should be a named vector indicating which samples belong to the  
# different groups.
doBartlett = function(inpMat, groups=NULL) {
  bartlett.p = apply(inpMat, 1, function(x) { tryCatch(bartlett.test(x, g=groups)$p.value, error = function(e) NA) })
  
  return(bartlett.p)
}

# Perform Levene's test on each row of the given matrix.
# groups is a factor indicating which samples belong to the different groups.
# location one of "median", "mean", "trim.mean"
doLevene = function(inpMat, groups, location=c("median", "mean", "trim.mean")) {
	# Load the package for performing Levene's test
	require(lawstat)
	# Get the right 
	location = match.arg(location)
	
	levene.p = apply(inpMat, 1, function(x) { tryCatch(levene.test(x, g=groups)$p.value, error = function(e) NA) })
  
	return(levene.p)
}

# Calculate prioritization score by summing over the rows
# x should be a matrix with data types in the columns
dataTypeScore = function(mat) {
  apply(mat, 1, sum, na.rm=TRUE) 
}

greaterThan = function(x, y) { x > y }
smallerThan = function(x, y) { x < y }

##' Creates a helper function to count the elements greater
##' than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the number of elements greater than the cutoff
##' @export
##' @author Andreas Schlicker
gtCutoff = function(cutoff) { function(x) { sum(x > cutoff) } }

##' Creates a helper function to find the percentage of elements
##' that are greater than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the percentage of elements greater than the cutoff (ranging from
##' 0 to 1.
##' @export
##' @author Andreas Schlicker
gtCutoffPercent = function(cutoff) { function(x) { sum(x > cutoff) / length(x) } }

##' Creates a helper function to count the elements less
##' than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the number of elements less than the cutoff
##' @export
##' @author Andreas Schlicker
ltCutoff = function(cutoff) { function(x) { sum(x < cutoff) } }

##' Creates a helper function to find the percentage of elements
##' that are less than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the percentage of elements less than the cutoff (ranging from
##' 0 to 1.
##' @export
##' @author Andreas Schlicker
ltCutoffPercent = function(cutoff) { function(x) { sum(x < cutoff) / length(x) } }

##' Convert a vector into a matrix with one row.
##' Names of the vector are preserved as colnames.
##' @param vec input vector
##' @return the new matrix
##' @author Andreas Schlicker
matrixFromVector = function(vec) {
	tmp = matrix(vec, nrow=1)
	colnames(tmp) = names(vec)
	tmp
}

##' Generates a boxplot from the given data.
##' @param group1 values for sample group1
##' @param group2 values for sample group2
##' @param group1.lab label for group1
##' @param group2.lab label for group2
##' @param xlabel x-axis label 
##' @param ylabel y-axis label
##' @param main main header for the plot
##' @param pvalue pvalue that will be added to the plot as text
boxplot = function(group1, group2, group1.lab, group2.lab, xlabel, ylabel, main, pvalue) {
	if (!require(ggplot2)) {
		stop("Can't load required package \"ggplot2\"!")
	}
	
	temp = data.frame(values=c(group1, group2), 
					  group=c(rep(group1.lab, times=length(group1)),
							  rep(group2.lab, times=length(group2))))
	
    cbbPalette = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	
	ggplot(temp, aes(x=group, y=values, fill=group)) +
	geom_boxplot(outlier.size=0, aes(alpha=0.3)) +
	geom_jitter(size=3) +
	guides(fill=FALSE, alpha=FALSE) +
	scale_colour_manual(values=cbbPalette) +
	geom_text(data=NULL, x=1, y=max(temp$values), label=paste("p = ", pvalue, sep="")) +
	xlab(xlabel) +
	ylab(ylabel) +
	ggtitle(main) + 
	theme(plot.title=element_text(face='bold', size=16),
		  panel.background = element_rect(fill='grey95', colour="grey95"),
		  axis.text.x=element_text(face='bold', size=16),
		  axis.title.x=element_text(face='bold', size=16),
		  axis.text.y=element_text(face='bold', size=16),
		  axis.title.y=element_text(face='bold', size=16),
		  legend.text=element_text(face="bold", size=16),
	      legend.title=element_text(face="bold", size=16),
		  strip.text.x = element_text(face="bold", size=16),
		  strip.text.y = element_text(face="bold", size=16))
}

