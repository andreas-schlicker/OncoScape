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
				   xlabel=NULL, ylabel=NULL, main=NULL, pvalue=NULL,
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
		generateTheme()
	
	if (!is.null(xlabel)) {
		p = p + xlab(xlabel)
	}
	if (!is.null(ylabel)) {
		p = p + ylab(ylabel)
	}
	if (!is.null(main)) {
		p = p + ggtitle(main)
	}
	if (!is.null(pvalue) && !is.na(pvalue)) {
		xpos = 1
		if (max(subset(temp, group==lab.group2)$values, na.rm=TRUE) > max(subset(temp, group==lab.group1)$values, na.rm=TRUE)) {
			xpos = 2
		}
		p = p + geom_text(x=xpos, y=max(temp$values, na.rm=TRUE), label=paste("FDR = ", signif(pvalue, digits=2), sep=""))
	}
	
	p
}

##' Generates a barplot from Project Achilles data. 
##' @param scores vector of phenotype scores from Project Achilles
##' @param upper.threshold threshold above which phenotype scores are considered to be significant
##' @param lower.threshold threshold below which phenotype scores are considered to be significant
##' @param main main header for the plot 
##' @return the plot
##' @author Andreas Schlicker
barplot = function(scores, upper.threshold=NULL, lower.threshold=NULL, main=NULL) {
	if (!require(ggplot2)) {
		stop("Can't load required package \"ggplot2\"!")
	}
	
	if (!require(stringr)) {
		stop("Can't load required package \"stringr\"!")
	}
	
	plotting.df = data.frame(phenoscore=scores, 
							 cls=sapply(names(scores), function(x) { str_sub(x, start=1, end=str_locate(x, "_")[1, 1]-1) }))
	
	p = ggplot(plotting.df, aes(x=cls, y=phenoscore, fill="#767676")) +
		geom_bar(stat="identity", position=position_dodge()) +
		geom_hline(yintercept=0, linetype=1) +
		xlab("Cell lines") +
		ylab("Achilles phenotype score") +
		guides(fill=FALSE) +
		generateTheme()

	if (!is.null(main)) {
		p = p + ggtitle(main)
	}
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
		scale_colour_manual(name="", values=color.palette) +
		ylim(0,1) +
		xlab("Methylation probes") +
		ylab("Average beta value") + 
		generateTheme() + 
		theme(axis.text.x=element_text(face='bold', size=16, angle=90))
	
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
	
	# Copy number plot
	samp1 = colnames(acgh.group1)
	samp2 = colnames(acgh.group2)
	if (!is.null(samples)) {
		samp1 = intersect(samples, colnames(acgh.group1))
		samp2 = intersect(samples, colnames(acgh.group2))
	}
	cn.box = boxplot(acgh.group1[gene, samp1], acgh.group2[gene, samp2], 
					 lab.group1, lab.group2, 
					 xlabel=NULL, ylabel=paste(gene, "copy number"), main=NULL, pvalue=prior.details[gene, "cgh.bh"],
					 color.palette)
			 
	# Achilles plot
	achil = barplot(achilles, achilles.ut, achilles.lt, main=NULL)
	
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
