# Some functions to generate plots for specific genes
# 
# Author: schlandi
###############################################################################

##' Generates a boxplot from the given data.
##' @param group1 values for sample group1
##' @param group2 values for sample group2
##' @param lab.group1 label for group1
##' @param lab.group2 label for group2
##' @param xlabel x-axis label 
##' @param ylabel y-axis label
##' @param main main header for the plot
##' @param pvalue pvalue that will be added to the plot as text
##' @param color.palette character vector with colors for plotting
##' @return the boxplot
##' @author Andreas Schlicker
boxplot = function(group1, group2, 
				   lab.group1="group1", lab.group2="group2", 
				   xlabel="", ylabel="", main=NULL, pvalue=NULL,
				   color.palette=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) {
	if (!require(ggplot2)) {
		stop("Can't load required package \"ggplot2\"!")
	}
	
	temp = data.frame(values=c(group1, group2), 
			group=c(rep(lab.group1, times=length(group1)),
					rep(lab.group2, times=length(group2))))
	
	p = ggplot(temp, aes(x=group, y=values, fill=group)) +
		geom_boxplot(outlier.size=0, aes(alpha=0.3)) +
		geom_jitter(size=3) +
		guides(fill=FALSE, alpha=FALSE) +
		scale_fill_manual(values=color.palette) +
		geom_hline(yintercept=0, linetype=1) +
		xlab(xlabel) + 
		ylab(ylabel) + 
		generateTheme()
	
	if (!is.null(main)) {
		p = p + ggtitle(main)
	}
	if (!is.null(pvalue) && !is.na(pvalue)) {
		xpos = 1
		if (max(subset(temp, group==lab.group2)$values, na.rm=TRUE) > max(subset(temp, group==lab.group1)$values, na.rm=TRUE)) {
			xpos = 2
		}
		p = p + geom_text(x=xpos, y=max(temp$values, na.rm=TRUE), label=paste("FDR = ", signif(pvalue, digits=2), sep=""), size=8, fontface="bold")
	}
	
	p
}

##' Generates a barplot from Project Achilles data. 
##' @param scores vector of phenotype scores from Project Achilles
##' @param upper.threshold threshold above which phenotype scores are considered to be significant
##' @param lower.threshold threshold below which phenotype scores are considered to be significant
##' @param main main title for the plot 
##' @return the plot
##' @author Andreas Schlicker
achillesBarplot = function(scores, upper.threshold=NULL, lower.threshold=NULL, main="") {
	require(ggplot2) || stop("Can't load required package \"ggplot2\"!")
	require(stringr) || stop("Can't load required package \"stringr\"!")
	
	plotting.df = data.frame(phenoscore=scores, 
							 cls=sapply(names(scores), function(x) { str_sub(x, start=1, end=str_locate(x, "_")[1, 1]-1) }))
	
	p = barplot(plotting.df, facet=NULL, x="cls", y="phenoscore", stat="identity", 
				title=main, xlab="Cell lines", ylab="Achilles phenotype score", fill="#767676") +
		geom_hline(yintercept=0, linetype=1) +
		theme(axis.text.x=element_text(face='bold', size=25, angle=45, vjust=0.5))
					 
	if (!is.null(upper.threshold)) {
		p = p + geom_hline(yintercept=upper.threshold, linetype=2)
	}
	if (!is.null(lower.threshold)) {
		p = p + geom_hline(yintercept=lower.threshold, linetype=2)
	}
	
	p
}

##' Generates a scatterplot from methylation data.
##' @param meth.group1 matrix of beta values, probes in columns, samples in rows
##' @param meth.group2 matrix of beta values, probes in columns, samples in rows
##' @param error.bar character string indicating statistic that is shown as error bars.
##' "se" for standard error, "ci" for confidence interval, "none" to prevent
##' plotting of error bars.
##' @param conf.interval percentage range for confidence interval; default: 0.95. This
##' parameter is only used if error bars represent confidence intervals.
##' @param lab.group1 label for group1
##' @param lab.group2 label for group2
##' @param main main header for the plot
##' @param color.palette colors for plotting
##' @return the scatterplot
##' @author Andreas Schlicker
scatterplot = function(meth.group1, meth.group2, 
					   error.bar=c("se", "ci", "none"), conf.interval=0.95,
					   lab.group1="group1", lab.group2="group2", main=NULL,
					   color.palette=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) {
				   
	error.bar = match.arg(error.bar)
	
	# Generate the combined data.frame
	data.df = data.frame(beta=c(as.vector(meth.group1), as.vector(meth.group2)),
					     probe=factor(c(rep(rownames(meth.group1), times=ncol(meth.group1)),
								        rep(rownames(meth.group2), times=ncol(meth.group2))),
									  levels=rownames(meth.group1)),
						 group=c(rep(lab.group1, times=nrow(meth.group1)*ncol(meth.group1)),
								 rep(lab.group2, times=nrow(meth.group2)*ncol(meth.group2))))

	# and summarize it for plotting
	plotting.df = summarySE(data.df, measurevar="beta", groupvars=c("probe", "group"),
							na.rm=TRUE, conf.interval=conf.interval)
	
	# Get the dodge to avoid error bar overplotting
	pd = position_dodge(0.5)
	
	p = ggplot(plotting.df, aes(x=probe, y=beta, shape=group, color=group)) +
		geom_point(size=4, position=pd) +
		scale_shape_manual(name="", values=c(15, 16)) +
		scale_color_manual(name="", values=color.palette) +
		ylim(0,1) +
		xlab("Methylation probes") +
		ylab("Average beta value") + 
		generateTheme() + 
		theme(axis.text.x=element_text(face='bold', size=25, angle=90))
	
	if (error.bar == "se") {
		p = p + geom_errorbar(data=plotting.df, aes(ymin=beta-se, ymax=beta+se), width=0.2, position=pd)
	} else if (error.bar == "ci") {
		p = p + geom_errorbar(data=plotting.df, aes(ymin=beta-ci, ymax=beta+ci), width=0.2, position=pd)
	}
	
	if (!is.null(main)) {
		p = p + ggtitle(main)
	}
	
	p
}

##' Creates all plots for one gene.
##' @param gene gene symbol
##' @param prior.details data.frame with detailed scores from the prioritization
##' @param samples optional character vector giving sample names to be used for
##' plotting
##' @param exprs.group1 expression matrix for group 1
##' @param exprs.group2 expression matrix for group 2
##' @param meth.group1 methylation matrix for group 1
##' @param meth.group2 methylation matrix for group 2
##' @param meth.anno annotation of methylation probes
##' @param acgh.group1 copy number matrix for group 1
##' @param acgh.group2 copy number matrix for group 2
##' @param achilles phenotype scores from project achilles
##' @param achilles.ut upper threshold used for scoring achilles
##' @param achilles.lt lower threshold used for scoring achilles
##' @param lab.group1 label for group 1; default: Tumors
##' @param lab.group2 label for group 2; default: Normals
##' @param color.palette colors to use for plotting
##' @return the combined plots
##' @author Andreas Schlicker
plotGene = function(gene, prior.details, samples=NULL, 
				 	exprs.group1, exprs.group2, 
		 			meth.group1, meth.group2, meth.anno, 
					acgh.group1, acgh.group2, 
					achilles, achilles.ut, achilles.lt, 
					lab.group1="Tumors", lab.group2="Normals", 
					color.palette=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) {

	# Gene expression plot
	samp1 = colnames(exprs.group1)
	if (!is.null(samples)) {
		samp1 = intersect(samples, colnames(exprs.group1))
	}
	ge.box = boxplot(exprs.group1[gene, intersect(samp1, intersect(colnames(exprs.group1), colnames(exprs.group2)))], 
					 exprs.group2[gene, intersect(samp1, intersect(colnames(exprs.group1), colnames(exprs.group2)))], 
					 lab.group1, lab.group2, 
					 xlabel=NULL, ylabel=paste(gene, "expression"), main=NULL, pvalue=prior.details[gene, "exprs.bh"],
					 color.palette=color.palette)
#	ge.box = boxplot(exprs.group1[gene, intersect(samp1, colnames(exprs.group1))], 
#			exprs.group2[gene, intersect(samp1, colnames(exprs.group2))], 
#			lab.group1, lab.group2, 
#			xlabel=NULL, ylabel=paste(gene, "expression"), main=NULL, pvalue=prior.details[gene, "exprs.bh"],
#			color.palette=color.palette)
	
	# Copy number plot
	samp1 = colnames(acgh.group1)
	samp2 = colnames(acgh.group2)
	if (!is.null(samples)) {
		samp1 = intersect(samples, colnames(acgh.group1))
		samp2 = intersect(samples, colnames(acgh.group2))
	}
	if (length(intersect(colnames(prior.details), "cgh.bh")) == 1) {
		pvalue = prior.details[gene, "cgh.bh"]
	} else { 
		pvalue = NA
	}
	cn.box = boxplot(acgh.group1[gene, samp1], acgh.group2[gene, samp2], 
					 lab.group1, lab.group2, 
					 xlabel=NULL, ylabel=paste(gene, "copy number"), main=NULL, pvalue=pvalue,
					 color.palette)
			 
	# Achilles plot
	achil = achillesBarplot(achilles, achilles.ut, achilles.lt, main="")
	
	# Methylation plot
	samp1 = colnames(meth.group1)
	samp2 = colnames(meth.group2)
	if (!is.null(samples)) {
		samp1 = intersect(samples, colnames(meth.group1))
		samp2 = intersect(samples, colnames(meth.group2))
	}
	meth.probes = rownames(meth.anno[which(meth.anno[, "genesym"] == gene), ])
	meth = scatterplot(meth.group1[meth.probes, samp1], meth.group2[meth.probes, samp2], 
					   error.bar="se", lab.group1=lab.group1, lab.group2=lab.group2, main=NULL,
					   color.palette=color.palette)
			   
	list(gene.expression=ge.box, acgh=cn.box, achilles=achil, methylation=meth)
}

##' Generate the standard theme for all plots.
##' @return the theme
##' @author Andreas Schlicker
generateTheme = function() {
	theme(plot.title=element_text(face='bold', size=25),
		  panel.background = element_rect(fill='grey95', color="grey95"),
		  axis.text.x=element_text(face='bold', size=25),
		  axis.title.x=element_text(face='bold', size=25),
		  axis.text.y=element_text(face='bold', size=25),
		  axis.title.y=element_text(face='bold', size=25),
		  legend.text=element_text(face="bold", size=25),
		  legend.title=element_text(face="bold", size=25),
		  strip.text.x = element_text(face="bold", size=25),
		  strip.text.y = element_text(face="bold", size=25))
}

## Taken from: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_%28ggplot2%29/
##
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
		conf.interval=.95, .drop=TRUE) {
	require(plyr)
	
	# New version of length which can handle NA's: if na.rm==T, don't count them
	length2 <- function (x, na.rm=FALSE) {
		if (na.rm) sum(!is.na(x))
		else       length(x)
	}
	
	# This is does the summary; it's not easy to understand...
	datac <- ddply(data, groupvars, .drop=.drop,
			.fun= function(xx, col, na.rm) {
				c( N    = length2(xx[,col], na.rm=na.rm),
						mean = mean   (xx[,col], na.rm=na.rm),
						sd   = sd     (xx[,col], na.rm=na.rm)
				)
			},
			measurevar,
			na.rm
	)
	
	# Rename the "mean" column   
	datac <- rename(datac, c("mean"=measurevar))
	
	datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
	
	# Confidence interval multiplier for standard error
	# Calculate t-statistic for confidence interval:
	# e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
	ciMult <- qt(conf.interval/2 + .5, datac$N-1)
	datac$ci <- datac$se * ciMult
	
	return(datac)
}

##' Sorts elements in the input data frame according to the given score in ascending order.
##' Used for plotting the results of a prioritization run.
##' @param inpData input data frame
##' @param score score type to use for sorting
##' @param column column name or index containing the elements to sort
##' @return character vector with sorted elements
##' @author Andreas Schlicker
getOrder = function(inpData, score, column) {
	# Retain only the sorting scores
	inpData = inpData[which(inpData[, "score.type"] == score), ]
	# Split the scores according to the different entries in "column"
	scores = split(inpData$score, inpData[, column])
	# Sum the scores
	res = unlist(lapply(scores, sum, na.rm=TRUE))
	# Return the elements sorted according to ascending sum of scores
	names(sort(res))
}

##' Summarize the results in the list for plotting as score heatmap.
##' @param results named list with prioritization results
##' @param scores named list of scores that are summarized;
##' The list elements have to match to column names in the results matrices.
##' The names of the list elements will be used as labels in the resulting data.frame.
##' Combined.affected is a special name, which is translated to og.affected.rel - ts.affected.rel
##' @return data.frame for plotting
##' @author Andreas Schlicker
heatmapDataframe = function(results, 
							scores=list(OG="og.score", TS="ts.score", Combined="combined.score", 
										OG.affected="og.affected.rel", TS.affected="ts.affected.rel", Combined.affected="")) {
	result.df = data.frame()
	for (n in names(results)) {
		temp = results[[n]]$prioritize.combined
		for (s in names(scores)) {
			if (s == "Combined.affected") {
				score = temp[, "og.affected.rel"] - temp[, "ts.affected.rel"]
			} else {
				score = temp[, scores[[s]]]	
			}
			result.df = rbind(result.df,
				data.frame(gene=rownames(temp),
						   score=score,
						   score.type=s,
						   cancer=n))
		}
	}
	
	result.df
}

##' Generate a heatmap of prioritization scores. Used for result heatmaps in the plot_prioritize script.
##' @param dataFrame approprately formatted data.frame
##' @param yaxis.theme ggplot2 theme for the y-axis
##' @param labels vector with labels for elements in rows; default: NULL (no labels)
##' @param breaks vector with y-axis breaks; is only used if labels != NULL; default: NULL
##' @param color.low color for the lowest value; default: white
##' @param color.mid color for middle values; default: NULL
##' @param color.high color for high values; default: black
##' @param title main title for the plot; default: ""
##' @param ylab title for the y-axis of the plot; default: ""
##' @param xlab title for the x-axis of the plot; default: ""
##' @return a ggplot2 object
##' @author Andreas Schlicker
getHeatmap = function(dataFrame, yaxis.theme, labels=NULL, breaks=NULL, color.low="white", color.mid=NULL, color.high="black", title="", ylab="", xlab="") {
	p = ggplot(dataFrame, aes(x=cancer, y=gene)) + 
			geom_tile(aes(fill=score), color = "white") + 
			labs(title=title, x=xlab, y=ylab) +
			theme(panel.background=element_rect(color="white", fill="white"),
					axis.ticks=element_blank(), 
					axis.text.x=element_text(color="grey50", face="bold")) +
			yaxis.theme
	
	if (!is.null(color.mid)) {
		p = p + scale_fill_gradient2(low=color.low, mid=color.mid, high=color.high)
	} else {
		p = p + scale_fill_gradient(low=color.low, high=color.high)
	}
	if (!is.null(labels)) {
		p = p + scale_y_discrete(breaks=breaks, labels=labels)
	}
	
	p
}

##' Generate a distribution plot for prioritization results. 
##' @param dataFrame approprately formatted data.frame
##' @param facets string that defines a formula for generating the facets of the plot
##' @param plot.type either "histogram" or "density"; default: "histogram"
##' @param ncol number of facets to plot per row; default: 3
##' @param title main title for the plot; default: ""
##' @param ylab title for the y-axis of the plot; default: ""
##' @param xlab title for the x-axis of the plot; default: ""
##' @param strip.theme ggplot2 theme for the facet strip; default: NULL
##' @return a ggplot2 object
##' @author Andreas Schlicker
getDistPlot = function(dataFrame, facets, plot.type=c("histogram", "density"), 
					   ncol=3, title="", xlab="", ylab="",
					   strip.theme=NULL) {
	plot.type = match.arg(plot.type)
	
	p = ggplot(dataFrame, aes(x=score, fill=score.type, color=score.type)) + 
			scale_fill_manual("Score", breaks=c("OG", "TS", "Combined"), values=c("#b52c2c", "#1575c6", "grey60")) +
			scale_color_manual("Score", breaks=c("OG", "TS", "Combined"), values=c("#b52c2c", "#1575c6", "grey60")) +
			scale_x_continuous(breaks=-4:4) + 
			facet_wrap(formula(facets), ncol=ncol) +
			labs(title=title, x=xlab, y=ylab) + 
			theme(axis.text=element_text(color="grey50", face="bold"))
	
	if (plot.type == "histogram") {
		p = p + geom_histogram(binwidth=0.5, position="dodge")
	} else {
		p = p + geom_density(size=1, alpha=0.3)
	}
	
	if (!is.null(strip.theme)) {
		p = p + strip.theme
	}
	
	p
}

##' Summarize the given score from the results in a matrix for plotting.
##' @param results the results list from the prioritization
##' @param score which score should be plotted; has to be a valid index or column name
##' for the prioritize.combined matrix; default: "combined.score"
##' @param summarize function for summarizing scores for each gene across cancers; defaul: NULL
##' @return the score matrix with genes in rows and cancers in columns
##' @author Andreas Schlicker
summaryMatrix = function(results, score, summarize=NULL) {
	plotting = matrix(NA, nrow=nrow(results[[1]]$prioritize.combined), ncol=length(results))
	colnames(plotting) = names(results)
	rownames(plotting) = rownames(results[[1]]$prioritize.combined)
	for (cancType in colnames(plotting)) {
		plotting[, cancType] = results[[cancType]]$prioritize.combined[rownames(plotting), score]
	}
	if (!is.null(summarize)) {
		plotting = cbind(plotting, apply(plotting, 1, summarize))
		colnames(plotting)[ncol(plotting)] = as.character(substitute(summarize))
	}
	
	plotting
}

##' Generate a boxplot showing specific scores for all genes across all cancer types.
##' The plot contains one facet per cancer type and one for the mean value. The genes
##' can be grouped to contrast different gene classes with each other.
##' @param results the results list from the prioritization
##' @param groups named list of vectors with gene IDs grouping genes together; default: NULL (no groups)
##' @param score which score should be plotted; has to be a valid index or column name
##' for the prioritize.combined matrix; default: "combined.score"
##' @param title main title for the plot; default: ""
##' @param ylab title for the y-axis of the plot; default: ""
##' @param xlab title for the x-axis of the plot; default: ""
##' @return a ggplot2 object
##' @author Andreas Schlicker
scoreBoxplot = function(results, groups=NULL, score="combined.score", title="", xlab="", ylab=score) {
	require(ggplot2) || stop("Couldn't load required package \"ggplot2\"!")
	require(grid) || stop("Couldn't load required package \"grid\"!")
	
	plotting = summaryMatrix(results, score, mean)
	
	if (is.null(groups)) {
		groups = list(all=rownames(results[[1]]$prioritize.combined))
	}
	
	plotting.df = data.frame()
	for (cancType in colnames(plotting)) {
		for (group in names(groups)) {
			genes = intersect(groups[[group]], rownames(results[[1]]$prioritize.combined))
			plotting.df = rbind(plotting.df,
					data.frame(Cancer=rep(cancType, times=length(genes)),
							Score=plotting[genes, cancType],
							Class=rep(group, times=length(genes))))
		}
	}	
	plotting.df[, 2] = as.numeric(plotting.df[, 2])

	p = ggplot(plotting.df, aes(x=Class, y=Score, fill=Class)) +
			geom_boxplot(outlier.size=0) +
			geom_jitter(size=1) +
			facet_wrap(~ Cancer, ncol=4) +
			labs(title=title, xlab=xlab, ylab=ylab) +
			theme(axis.title=element_text(color="grey50", face="bold", size=20),
					axis.text=element_text(color="grey50", face="bold", size=20),
					legend.text=element_text(color="grey50", face="bold", size=20),
					legend.title=element_text(color="black", face="bold", size=20),
					legend.key.size=unit(1, "cm"),
					strip.text=element_text(color="grey20", face="bold", size=20))
	
	p
}

##' Generate a histogram visualizing specific scores for all genes belonging to certain groups.
##' The plot contains one facet per cancer type and one for the sum. The genes
##' can be grouped to contrast different gene classes with each other and genes that don't belong to any group.
##' At least one group is required.
##' @param results the results list from the prioritization
##' @param groups named list of vectors with gene IDs grouping genes together
##' @param score which score should be plotted; has to be a valid index or column name
##' for the prioritize.combined matrix; default: "combined.score"
##' @param title main title for the plot; default: ""
##' @param ylab title for the y-axis of the plot; default: ""
##' @param xlab title for the x-axis of the plot; default: ""
##' @param cols a color palette
##' @return a ggplot2 object
##' @author Andreas Schlicker
scoreHistogram = function(results, groups, score="combined.score", title="", xlab="", ylab=score, 
						  cols=c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) { 
	require(ggplot2) || stop("Couldn't load required package \"ggplot2\"!")
	
	plotting = summaryMatrix(results, score, sum)
	
	unique.genes = setdiff(rownames(plotting), unique(unlist(groups)))
	
	somscore.df = data.frame()
	for (cancType in colnames(plotting)) {
		for (group in names(groups)) {
			cgst = table(plotting[groups[[group]], cancType])
			ugst = table(plotting[unique.genes, cancType])
			
			somscore.df = rbind(somscore.df,
								data.frame(Cancer=rep(cancType, times=length(cgst)),
										   Score=names(cgst),
										   Percentage=round(cgst/length(common.genes), digits=2),
										   Group=rep(group, times=length(cgst))),
								data.frame(Cancer=rep(cancType, times=length(ugst)),
										   Score=names(ugst),
										   Percentage=round(ugst/length(unique.genes), digits=2),
										   Group=rep(group, times=length(ugst))))
		}
	}
	colnames(somscore.df)[3] = "Percentage"
	somscore.df[, 3] = as.numeric(somscore.df[, 3])
	somscore.df[, 2] = as.numeric(somscore.df[, 2])
	
	p = ggplot(somscore.df, aes(x=Score, y=Percentage, fill=Group)) + 
		geom_histogram(binwidth=0.1, position="dodge", stat="identity") + 
		facet_wrap( ~ Cancer, ncol=4) +
		scale_color_manual(values=cols) +
		labs(title=title, ylab=ylab, xlab=xlab) + 
		theme(axis.ticks=element_blank(), 
			  axis.text.x=element_text(color="grey50"),
			  axis.text.y=element_text(color="grey50"))
	
	p
}

##' Plots a confusion heatmap. Essentially a colored confusion table. The plot can be faceted.
##' @param dataframe plotting dataframe with at least four columns (x values, y values, frequencies and facet variable)
##' @param facet name of column to be used for faceting; if == NULL, no faceting is applied; default: column four
##' @param ncol number of facets in each row; default: 3
##' @param title main title for the plot; default: name of the third column of dataframe
##' @param xlab x-axis title; default: name of the first column of dataframe
##' @param ylab y-axis title; default: name of the second column of dataframe
##' @author Andreas Schlicker (adopted from function ggfluctuation() in ggplot2
confusionHeatmap = function(dataframe, facet=colnames(dataframe)[4], ncol=3,
							title=colnames(dataframe)[3], xlab=colnames(dataframe)[1], ylab=colnames(dataframe)[2]) {
	require(ggplot2) || stop("Can't load required package \"ggplot2\"!")
	require(stringr) || stop("Can't load required package \"stringr\"!")
						
	names(dataframe)[1:3] = c("x", "y", "freq")
	dataframe[, "x"] = factor(dataframe[, "x"], levels=sort(as.integer(levels(dataframe[, "x"]))))
	dataframe[, "y"] = factor(dataframe[, "y"], levels=sort(as.integer(levels(dataframe[, "y"]))))
	
	if (is.null(facet)) {
		temp = unlist(lapply(split(dataframe$freq, list(dataframe$x, dataframe$y)), sum))
		dataframe = dataframe[1:length(temp), ]
		for (i in 1:length(temp)) {
			dataframe[i, 3] = temp[i]
			dataframe[i, 1:2] = str_split(names(temp)[i], "\\.")[[1]]
		}
	}
	
	p = ggplot(dataframe, aes(x=x, y=y, fill=freq, label=freq)) +
		geom_tile(color = "grey50") +
		geom_text(size=10, fontface="bold", color="gray10") + 
		scale_fill_gradient2(name="Frequency", low="white", high="#35a435") +
		labs(title=title, x=xlab, y=ylab) + 
		theme(title=element_text(color="black", size=20, face="bold"),
			  axis.text=element_text(color="grey50", size=20, face="bold"),
			  axis.title=element_text(color="grey50", size=20, face="bold"),
			  legend.text=element_text(color="grey50", size=20, face="bold"),
			  legend.title=element_text(color="grey50", size=20, face="bold"),
			  strip.text=element_text(color="black", size=20, face="bold"))

	if (!is.null(facet)) {
		p = p + facet_wrap(formula(paste(" ~ ", facet, sep="")), ncol=ncol)
	}
	
	p
}

##' Creates a barplot for a specific column in the input data. 
##' By default, a facet will be created for each level of the "cancer" column.
##' @param dataframe input data
##' @param facet name of the column used for faceting; default: cancer
##' @param ncol number of columns for faceting; default: 3
##' @param x name of the column to be used as x-axis; default: score.diff
##' @param y name of the column to be used as y-axis; default: NULL (for use with stat=="bin")
##' @param stat statistic to be used for the plot; default: bin
##' @param title main title for the plot; default: ""
##' @param xlab x-axis title; default: ""
##' @param ylab y-axis title; default: ""
##' @param fill bar color; default: gray30 
##' @author Andreas Schlicker
barplot = function(dataframe, facet="cancer", ncol=3, 
				   x="score.diff", y=NULL, stat="bin",			   
				   title="", xlab="", ylab="", fill="gray30") {
    require(ggplot2) || stop("Can't load required package \"ggplot2\"!")
	
	dataframe[, x] = as.factor(dataframe[, x])
	
	if (is.null(y)) {
		p = ggplot(dataframe, aes_string(x=x), fill=fill)
	} else {
		p = ggplot(dataframe, aes_string(x=x, y=x), fill=fill)
	}
	
	p = p + geom_bar(binwidth=0.25, stat=stat) + 
			guides(fill=FALSE) +
			labs(title=title, xlab=xlab, ylab=ylab) +
			theme(title=element_text(color="black", size=20, face="bold"),
				  axis.text=element_text(color="grey50", size=20, face="bold"),
			  	  axis.title=element_text(color="grey50", size=20, face="bold"),
			  	  legend.text=element_text(color="grey50", size=20, face="bold"),
			  	  legend.title=element_text(color="grey50", size=20, face="bold"),
			  	  strip.text=element_text(color="black", size=20, face="bold"))

	if (!is.null(facet)) {
		p = p + facet_wrap(formula(paste(" ~ ", facet, sep="")), ncol=ncol)
	}

	p
}

##' Creates a data.frame that is suitable for plotting with confusionHeatmap.
##' @param dataframe a data.frame that contains the scores to be compared (column names
##' score and score2) and a cancer column
##' @return the new data.frame
##' @author Andreas Schlicker
conHeatDf = function(dataframe) {
	res.df = data.frame()
	for (cancType in unique(dataframe[, "cancer"])) {
		tempTab = table(subset(dataframe, cancer == cancType)[, c("score", "score2")])
		res.df = rbind(res.df, data.frame(x=rep(rownames(tempTab), times=ncol(tempTab)),
										  y=rep(colnames(tempTab), each=nrow(tempTab)),
										  freq=as.vector(tempTab),
										  cancer=cancType))
	}
	
	res.df
}

##' Multiple plot function
##' Taken from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_%28ggplot2%29/
##'
##' ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
##' @param cols Number of columns in layout
##' @param layout A matrix specifying the layout. If present, 'cols' is ignored.
##'
##' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
##' then plot 1 will go in the upper left, 2 will go in the upper right, and
##' 3 will go all the way across the bottom.
##' @author Winston Chang
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
	require(grid)
	
	# Make a list from the ... arguments and plotlist
	plots <- c(list(...), plotlist)
	
	numPlots = length(plots)
	
	# If layout is NULL, then use 'cols' to determine layout
	if (is.null(layout)) {
		# Make the panel
		# ncol: Number of columns of plots
		# nrow: Number of rows needed, calculated from # of cols
		layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
				ncol = cols, nrow = ceiling(numPlots/cols))
	}
	
	if (numPlots==1) {
		print(plots[[1]])
		
	} else {
		# Set up the page
		grid.newpage()
		pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
		
		# Make each plot, in the correct location
		for (i in 1:numPlots) {
			# Get the i,j matrix positions of the regions that contain this subplot
			matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
			
			print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
							layout.pos.col = matchidx$col))
		}
	}
}

##' Prints a heatmap that is defined by the given list of parameters.
##' @param params named list with parameters: "data": data frame for plotting; 
##' "cancer.order": sorted list of column names; "gene.order": sorted list of row names;
##' "topgenes": rowa to plot; "prefix": filename prefix including directory
##' "filename": filename including extension; "score.type": score to plot; "color.low": color for low values
##' "color.mid": color for intermediate values; "color.high": color for high values; 
##' "yaxis": ggplot2 theme to apply to the yaxis; "ylab": y-axis title; "title": plot title
##' @return NULL
##' @author Andreas Schlicker
plotHeatmap = function(params) {
	params$data[, "cancer"] = factor(params$data[, "cancer"], levels=params$cancer.order)
	params$data[, "gene"] = factor(params$data[, "gene"], levels=params$gene.order)
	
	png(paste(configOptions$prefix, params$filename, sep="_"), width=5000, height=3000, res=300)
	print(getHeatmap(dataFrame=subset(params$data, score.type==params$score.type & gene %in% params$topgenes), 
					yaxis.theme=params$yaxis, 
					color.low=params$color.low, color.mid=params$color.mid, color.high=params$color.high, 
					ylab=params$ylab, title=params$title))
	invisible(dev.off())
}

##' Summarizes the given score in a matrix. 
##' @param results the result list from a prioritization run
##' @param score name of the column to add to the matrix; default: combined.score
##' @return matrix with genes in rows and cancers in columns
##' @author Andreas Schlicker
pathviewMat = function(results, score="combined.score") {
	allGenes = c()
	for (i in 1:length(results)) {
		allGenes = union(allGenes, rownames(results[[i]]$prioritize.combined))
	}
	
	res = matrix(NA, nrow=length(allGenes), ncol=length(results))
	rownames(res) = allGenes
	colnames(res) = names(results)
	for (n in names(results)) {
		temp = results[[n]]$prioritize.combined[, score, drop=FALSE]
		names(temp) = rownames(results[[n]]$prioritize.combined)
		res[, n] = temp[allGenes]
	}
	
	res
}

##' Generates pathview plots for all given results. This function can be used to plot 
##' different scores for one cancer type or one score across different cancer types.
##' If cancers == "all", only the first score in scores is plotted. In all other cases,
##' all scores are plotted for the first cancer only. 
##' @param results the result list from a prioritization run
##' @param pathways vector with KEGG pathway IDs to plot; default: NULL (all pathways)
##' @param cancers vector of names of cancer types to plot; default: all
##' @param scores vector of names of the scores to plot; default: combined.score
##' @param combine function for combining scores across all cancer for each gene (e.g. min, max);
##' the function has to accept a logical na.rm parameter; default: NULL (do not combine, one plot per cancer)
##' @param out.dir directory for output files; default: "."
##' @param out.suffix suffix to be added to output plots; default: ""
##' @param kegg.dir directory with predownloaded KEGG files; all new downloaded files will be stored there;
##' @param multi.state logical, combine different samples or scores into one plot; default: FALSE
##' default: "."
##' @param multiplier factor to multiply all scores with; can be used to transform tumor suppressor genes to 
##' negative values; default: 1
##' @param colors named list of colors for plotting gene values; if specified, has to contain elements for "low",
##' "mid" and "high"; default: NULL
##' @author Andreas Schlicker
generatePathview = function(results, pathways, cancers="all", scores="combined.score", combine=NULL,
							out.dir=".", out.suffix="", kegg.dir=".", multi.state=FALSE, multiplier=1,
							colors=NULL) {
	require(pathview) || stop("Can't load required package \"pathview\"!")
	
	# Get the score matrix and combine them if necessary
	if (cancers == "all") {
		cancers = names(results)
	}
	
	if (length(cancers) > 1) {
		scoreMat = pathviewMat(results[cancers], scores[1])
	} else {
		scoreMat = results[[cancers[1]]]$prioritize.combined[, scores]
	}
	if (!is.null(combine)) {
		scoreMat = apply(scoreMat, 1, combine, na.rm=TRUE)
	}
	scoreMat = scoreMat * multiplier
	
	# Save working directory and switch to new one
	oldwd = getwd()
	setwd(out.dir)
	
	low = list(gene="green", cpd="blue")
	mid = list(gene="gray", cpd="gray")
	high = list(gene="red", cpd="yellow")
	if (!is.null(colors)) {
		low$gene = colors$low
		mid$gene = colors$mid
		high$gene = colors$high
	}
	
	# Generate plots
	for (path in pathways) {
		invisible(pathview(gene.data=scoreMat, 
				 		   pathway.id=path,
				 		   kegg.native=TRUE, 
		   				   gene.idtype="SYMBOL", 
			   	 		   limit=list(gene=max(abs(scoreMat), na.rm=TRUE), cpd=1), 
			 		   	   node.sum="max.abs", 
			 		   	   out.suffix=out.suffix,
			 		   	   kegg.dir=kegg.dir,
						   multi.state=multi.state,
						   low=low, mid=mid, high=high))
   	}
	 
	# Back to old directory
	setwd(oldwd)
}
