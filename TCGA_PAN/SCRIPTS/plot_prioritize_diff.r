#!/usr/bin/Rscript

library(ggplot2)
library(grid)

source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/utils.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/compare.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/plotting.r")

OPTIONS = getOptionList(c("config", "prefix", "result1", "result2"))

# Is there a configuration file?
args = commandArgs(TRUE)
if (length(args) == 0) {
	args = c("--help")
}
config = parseOptions(OPTIONS, args)
if (!is.null(config$config)) {
	config = parseConfigFile(config$config)
} else {
	config = c()
}
# Parse options from configuration file and command line arguments
configOptions = parseOptions(OPTIONS, c(config, commandArgs(TRUE)))


# Load prioritization results 1
load(configOptions$result1)
results1 = results

# Load prioritization results 2
load(configOptions$result2)
results2 = results

for (n in names(results1)) {
	commonGenes = intersect(rownames(results1[[n]]$prioritize.combined), rownames(results2[[n]]$prioritize.combined))
	results1[[n]]$prioritize.combined = results1[[n]]$prioritize.combined[commonGenes, ]
	results2[[n]]$prioritize.combined = results2[[n]]$prioritize.combined[commonGenes, ]
}

results.diff = diffMat(results1, results2)

result.df = heatmapDataframe(results.diff)

parameters = list()

## TS heatmaps
gene.order = getOrder(result.df, "TS", "gene")
diffGenes = names(alteredGenes(result.df, "TS"))
diffGenes = c(diffGenes[1:50], diffGenes[(length(diffGenes)-50+1):length(diffGenes)])
cancer.order = getOrder(result.df, "TS", "cancer")
parameters[["ts3"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_heatmap.png", gene.order=gene.order,
													 topgenes=diffGenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="grey50", face="bold")), score.type="TS", 
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Difference in tumor suppressor score")
parameters[["ts4"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_affected_heatmap.png", gene.order=gene.order,
													 topgenes=diffGenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="grey50", face="bold")), score.type="TS.affected", 
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Difference in % affected samples TS")

## OG heatmaps
gene.order = getOrder(result.df, "OG", "gene")
diffGenes = names(alteredGenes(result.df, "OG"))
diffGenes = c(diffGenes[1:50], diffGenes[(length(diffGenes)-50+1):length(diffGenes)])
cancer.order=getOrder(result.df, "OG", "cancer")
parameters[["og3"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_heatmap_topgenes.png", gene.order=gene.order,
													 topgenes=diffGenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="grey50", face="bold")), score.type="OG",
													 color.low="#034b87", color.mid="gray98", color.high="#880000", 
													 ylab="", title="Difference in oncogene score")
parameters[["og4"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_affected_heatmap.png", gene.order=gene.order,
													 topgenes=diffGenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="grey50", face="bold")), score.type="OG.affected",
													 color.low="#034b87", color.mid="gray98", color.high="#880000", 
													 ylab="", title="Difference in % affected samples OG")


## Combined score heatmaps
gene.order = getOrder(result.df, "Combined", "gene")
diffGenes = names(alteredGenes(result.df, "Combined"))
diffGenes = c(diffGenes[1:50], diffGenes[(length(diffGenes)-50+1):length(diffGenes)])
cancer.order=getOrder(result.df, "Combined", "cancer")
parameters[["cs3"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_heatmap.png", gene.order=gene.order,
													 topgenes=diffGenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="grey50", face="bold")), score.type="Combined",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Difference in combined score")
parameters[["cs4"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_affected_heatmap.png", gene.order=gene.order,
													 topgenes=diffGenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="grey50", face="bold")), score.type="Combined.affected",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Difference in % affected samples Combined score")

# And plot the heatmaps		 
invisible(lapply(parameters, plotHeatmap))

## Score distributions
result.df[, "gene"] = factor(result.df[, "gene"], levels=gene.order)
result.df[, "cancer"] = factor(result.df[, "cancer"], levels=cancer.order)
png(paste(configOptions$prefix, "score_histogram_by_cancer.png", sep="_"), width=4000, height=3000, res=300)
print(getDistPlot(subset(result.df, score.type %in% c("OG", "TS", "Combined")), facets="~ cancer",
									plot.type="histogram", ncol=3,
									title="Score distribution by cancer", x="Score", y="Number of genes"))
invisible(dev.off())

png(paste(configOptions$prefix, "score_density_by_cancer.png", sep="_"), width=4000, height=3000, res=300)
print(getDistPlot(subset(result.df, score.type %in% c("OG", "TS", "Combined")), facets="~ cancer",
									plot.type="density", ncol=3,
									title="Score distribution by cancer", x="Score", y="Density"))
invisible(dev.off())



scores = c("Copy number", "Methylation", "Mutation", "Achilles", "Expression")
names(scores) = c("cna", "methylation", "mutations", "achilles", "exprs")
score.colors = c("CNA"="#999999", "Expr"="#E69F00", "Meth"="#56B4E9", "Mut"="#009E73", "shRNA"="#F0E442", "others"="#CC79A7")
score.colors2 = score.colors
names(score.colors2) = c("Copy number", "Expression", "Methylation", "Mutation", "Achilles", "others")


cancer.score.df = data.frame()
for (cancType in names(results.diff)) {
	for (score in names(scores)) {
		for (st in c("og", "ts")) {
			cancer.score.df = rbind(cancer.score.df,
															data.frame(Cancer=cancType,
																		 		 Gene=rownames(results.diff[[cancType]]$prioritize.combined),
																		 		 Category=scores[score],
																		 		 Score.category=toupper(st),
																		 		 Score=results.diff[[cancType]]$prioritize.combined[, paste(st, score, sep=".")]))
		}
	}
}

png(paste(configOptions$prefix, "category_distribution.png", sep="_"), width=3000, height=3000, res=300)
ggplot(cancer.score.df, aes(x=Score, fill=Category)) + 
geom_histogram(binwidth=0.4, position="dodge") + 
facet_grid(Cancer ~ Score.category) +
scale_x_continuous(breaks=-1:1) +
scale_fill_manual(values=score.colors2) + 
labs(title="Score distribution by gene", y="Number of genes", x="Score") + 
theme(axis.ticks=element_blank(), 
		 	axis.text.x=element_text(colour="grey50", face="bold"),
		 	axis.text.y=element_text(colour="grey50", face="bold"))
invisible(dev.off())


# Plot confusion heatmaps and barplots for single scores
for (m in c("combined.score", "og.score", "og.cna", "og.methylation", "og.mutations", "og.exprs", 
						"ts.score", "ts.cna", "ts.methylation", "ts.mutations", "ts.exprs")) {
	notEqual = sapply(names(results1), function(n) { any(results1[[n]]$prioritize.combined[, m] - 
																											 results2[[n]]$prioritize.combined[, m] != 0) })
	if (any(notEqual, na.rm=TRUE)) {
		plots = compareScore(results1, results2, list(OG=m))
		if (!(m %in% c("combined.score", "og.score", "ts.score"))) {
			png(paste(configOptions$prefix, "_comparison_confusion_", m, ".png", sep=""), width=4000, height=4000, res=300)
			print(plots[[1]])
			dev.off()
		}
		png(paste(configOptions$prefix, "_comparison_confusion_", m, "_sum.png", sep=""), width=4000, height=4000, res=300)
		print(plots[[3]])
		dev.off()
		png(paste(configOptions$prefix, "_comparison_barplot_", m, ".png", sep=""), width=4000, height=4000, res=300)
		print(plots[[2]])
		dev.off()
		png(paste(configOptions$prefix, "_comparison_barplot_", m, "_sum.png", sep=""), width=4000, height=4000, res=300)
		print(plots[[4]])
		dev.off()
	}
}

