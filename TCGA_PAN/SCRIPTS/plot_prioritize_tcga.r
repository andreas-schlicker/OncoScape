#!/usr/bin/Rscript

library(ggplot2)
library(grid)

source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/utils.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/plotting.r")

OPTIONS = getOptionList(c("config", "input", "prefix", "notopgenes"))

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


# Load prioritization results
load(configOptions$input)

MAXGENES = min(1000, nrow(results[[1]]$prioritize.combined))

result.df = heatmapDataframe(results)
result.df[which(result.df[, "score.type"] == "TS"), "score"] = -1 * result.df[which(result.df[, "score.type"] == "TS"), "score"]

## Fill the list of parameters
parameters = list()

gene.order = getOrder(result.df, "TS", "gene")
topgenes = gene.order[1:min(MAXGENES, length(gene.order))]
cancer.order = getOrder(result.df, "TS", "cancer")
parameters[["ts1"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_blank()), score.type="TS", 
													 color.low="#034b87", color.high="gray98",
													 ylab="", title="Tumor suppressor score")
parameters[["ts2"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_affected_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_blank()), score.type="TS.affected", 
													 color.low="gray98", color.high="#034b87",
													 ylab="", title="% affected samples TS")
													 
topgenes=gene.order[1:(2*configOptions$notopgenes)]
parameters[["ts3"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_heatmap_topgenes.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=20)), score.type="TS", 
													 color.low="#034b87", color.high="gray98",
													 ylab="", title="Tumor suppressor score")
parameters[["ts4"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_affected_heatmap_topgenes.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=20)), score.type="TS.affected", 
													 color.low="gray98", color.high="#034b87",
													 ylab="", title="% affected samples TS")
													 
gene.order = getOrder(result.df, "OG", "gene")
topgenes = gene.order[max(1, (length(gene.order)-MAXGENES-1)):length(gene.order)]
cancer.order = getOrder(result.df, "OG", "cancer")
parameters[["og1"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_blank()), score.type="OG",
													 color.low="gray98", color.high="#880000", 
													 ylab="", title="Oncogene score")
parameters[["og2"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_affected_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_blank()), score.type="OG.affected",
													 color.low="gray98", color.high="#880000", 
													 ylab="", title="% affected samples OG")

topgenes = gene.order[(length(gene.order)-(2*configOptions$notopgenes)+1):length(gene.order)]
parameters[["og3"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_heatmap_topgenes.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=20)), score.type="OG",
													 color.low="gray98", color.high="#880000", 
													 ylab="", title="Oncogene score")
parameters[["og4"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_affected_heatmap_topgenes.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=20)), score.type="OG.affected",
													 color.low="gray98", color.high="#880000", 
													 ylab="", title="% affected samples OG")
													 
gene.order=getOrder(result.df, "Combined", "gene")
topgenes= if (MAXGENES == length(gene.order)) {	
						topgenes = gene.order 
					} else { 
						topgenes = c(gene.order[1:ceiling(MAXGENES/2)], gene.order[(length(gene.order)-(floor(MAXGENES/2)-1)):length(gene.order)]) 
					} 
cancer.order = getOrder(result.df, "Combined", "cancer")
parameters[["cs1"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_blank()), score.type="Combined",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Combined score")
parameters[["cs2"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_affected_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_blank()), score.type="Combined.affected",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="% affected samples Combined score")
													 
topgenes = c(gene.order[1:configOptions$notopgenes], gene.order[(length(gene.order)-configOptions$notopgenes+1):length(gene.order)])
parameters[["cs3"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_heatmap_topgenes.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=20)), score.type="Combined",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Combined score")
parameters[["cs4"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_affected_heatmap_topgenes.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=20)), score.type="Combined.affected",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="% affected samples Combined score")
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

png(paste(configOptions$prefix, "score_histogram_by_gene.png", sep="_"), width=4000, height=3000, res=300)
print(getDistPlot(subset(result.df, score.type %in% c("OG", "TS", "Combined") & gene %in% topgenes), facets="~ gene",
									plot.type="histogram", ncol=4,
									title="Score distribution by gene", x="Score", y="Number of cancers"))
invisible(dev.off())

png(paste(configOptions$prefix, "score_density_by_gene.png", sep="_"), width=4000, height=3000, res=300)
print(getDistPlot(subset(result.df, score.type %in% c("OG", "TS", "Combined") & gene %in% topgenes), facets="~ gene",
									plot.type="density", ncol=4,
									title="Score distribution by gene", x="Score", y="Density"))
invisible(dev.off())


scores = c("Copy number", "Methylation", "Mutation", "Achilles", "Expression")
names(scores) = c("cna", "methylation", "mutations", "achilles", "exprs")
score.colors = c("CNA"="#999999", "Expr"="#E69F00", "Meth"="#56B4E9", "Mut"="#009E73", "shRNA"="#F0E442", "others"="#CC79A7")
score.colors2 = score.colors
names(score.colors2) = c("Copy number", "Expression", "Methylation", "Mutation", "Achilles", "others")

cancer.score.df = data.frame()
for (cancType in names(results)) {
	for (score in names(scores)) {
		for (st in c("og", "ts")) {
			cancer.score.df = rbind(cancer.score.df,
															data.frame(Cancer=cancType,
																		 		 Gene=rownames(results[[cancType]]$prioritize.combined),
																		 		 Category=scores[score],
																		 		 Score.category=ifelse(st=="og", "Oncogene", "Tumor Suppressor"),
																		 		 Score=results[[cancType]]$prioritize.combined[, paste(st, score, sep=".")]))
		}
	}
}


png(paste(configOptions$prefix, "category_distribution.png", sep="_"), width=3000, height=3000, res=300)
ggplot(subset(cancer.score.df, Score=="1"), aes(x=Score, fill=Category)) + 
geom_histogram(binwidth=0.4, position="dodge") + 
facet_grid(Cancer ~ Score.category) +
scale_x_continuous(breaks=0:1) +
scale_fill_manual(values=score.colors2) + 
labs(title="Score distribution by gene", y="Number of genes", x="Score") + 
theme(axis.ticks=element_blank(), 
		 	axis.text.x=element_text(colour="gray30", face="bold"),
		 	axis.text.y=element_text(colour="gray30", face="bold"))
invisible(dev.off())


## Summary of all cancer types
allGenes = unique(unlist(lapply(results, function(x) { rownames(x$prioritize.combined) })))

score.mat = matrix(0, nrow=length(allGenes), ncol=10)
colnames(score.mat) = paste(c("og.", "ts."), rep(c("cna", "exprs", "methylation", "mutations", "achilles"), each=2), sep="")
rownames(score.mat) = allGenes
score.cat = rep(c("CNA", "Expr", "Meth", "Mut", "shRNA"), each=2)
names(score.cat) = colnames(score.mat)

for (n in 1:length(results)) {
	for (score in colnames(score.mat)) {
		temp = names(which(results[[n]]$prioritize.combined[, score] == 1))
		score.mat[temp, score] = 1
	}
}

for (i in colnames(score.mat)) {
	cancer.score.df = rbind(cancer.score.df,
									 				data.frame(Cancer="any",
									 									 Gene=rownames(score.mat),
									 									 Category=score.cat[i],
									 									 Score.category=ifelse(length(grep("og.", i)) == 1, "Oncogene", "Tumor Suppressor"),
									 									 Score=score.mat[, i]))
}


cancer.score.df$Category = as.character(cancer.score.df$Category)
cancer.score.df[which(cancer.score.df$Category == "Copy number"), "Category"] = "CNA"
cancer.score.df[which(cancer.score.df$Category == "Expression"), "Category"] = "Expr"
cancer.score.df[which(cancer.score.df$Category == "Methylation"), "Category"] = "Meth"
cancer.score.df[which(cancer.score.df$Category == "Mutation"), "Category"] = "Mut"
cancer.score.df[which(cancer.score.df$Category == "Achilles"), "Category"] = "shRNA"
cancer.score.df$Category = factor(cancer.score.df$Category, level=sort(unique(as.character(cancer.score.df$Category))))

cancer.score.df$Score.category = factor(cancer.score.df$Score.category, level=c("Oncogene", "Tumor Suppressor"))
cancer.score.df$Cancer = factor(cancer.score.df$Cancer, level=sort(levels(cancer.score.df$Cancer)))
#cancer.score.df$Cancer = factor(cancer.score.df$Cancer, levels=c(levels(cancer.score.df$Cancer)[2:12], "any"))

png(paste(configOptions$prefix, "category_distribution2.png", sep="_"), width=4000, height=2500, res=300)
ggplot(subset(cancer.score.df, Score=="1"), aes(x=Category, fill=Score.category)) + 
geom_histogram(binwidth=0.3, position="dodge") + 
facet_wrap( ~ Cancer, nrow=4) +
scale_fill_manual(values=c("#D55E00", "#0072B2")) + 
labs(y="Number of genes", x="Score") + 
theme(text=element_text(color="gray30", size=20, face="bold"),
		 	title=element_text(color="gray30", size=20, face="bold"),
		 	axis.ticks=element_blank(), 
		 	legend.title=element_blank(),
		 	legend.key.size=unit(20, "points"))
invisible(dev.off())


png(paste(configOptions$prefix, "category_distribution_og.png", sep="_"), width=3500, height=2000, res=300)
ggplot(subset(cancer.score.df, Score=="1" & Score.category=="Oncogene"), aes(x=Category, fill=Category)) + 
geom_histogram(binwidth=0.3, position="dodge") + 
facet_wrap( ~ Cancer, ncol=4) +
#coord_flip() +
guides(fill=FALSE) + 
scale_fill_manual(values=score.colors) +
#scale_y_continuous(breaks=seq(2000, 6000, by=2000), labels=seq(2, 6, by=2)) +
labs(y="Number of genes", x="Score") + 
theme(text=element_text(color="gray30", size=20, face="bold"),
		 	title=element_text(color="gray30", size=20, face="bold"),
		 	legend.title=element_blank(),
		 	legend.key.size=unit(20, "points"))
invisible(dev.off())

png(paste(configOptions$prefix, "category_distribution_og_byscore.png", sep="_"), width=3500, height=1500, res=300)
ggplot(subset(cancer.score.df, Score=="1" & Score.category=="Oncogene"), aes(x=Cancer, fill=Category)) + 
geom_histogram(binwidth=0.3, position="dodge") + 
facet_wrap( ~ Category, ncol=5) +
coord_flip() +
guides(fill=FALSE) + 
scale_fill_manual(values=score.colors) +
scale_y_continuous(breaks=c(2000, 6000, 10000), labels=c(2000, 6000, 10000)) +
labs(y="Number of genes", x="") + 
theme(axis.text.x=element_text(color="gray30", size=16, face="bold", vjust=0),
			axis.text.y=element_text(color="gray30", size=20, face="bold"),
		 	title=element_text(color="gray30", size=20, face="bold"),
		 	strip.text=element_text(color="gray30", size=20, face="bold"),
		 	legend.title=element_blank())
dev.off()

png(paste(configOptions$prefix, "category_distribution_ts.png", sep="_"), width=3500, height=2000, res=300)
ggplot(subset(cancer.score.df, Score=="1" & Score.category=="Tumor Suppressor"), aes(x=Category, fill=Category)) + 
geom_histogram(binwidth=0.3, position="dodge") + 
facet_wrap( ~ Cancer, ncol=4) +
coord_flip() +
guides(fill=FALSE) + 
scale_fill_manual(values=score.colors) +
#scale_y_continuous(breaks=seq(2000, 6000, by=2000), labels=seq(2, 6, by=2)) +
labs(y="Number of genes", x="Score") + 
theme(text=element_text(color="gray30", size=20, face="bold"),
		 	title=element_text(color="gray30", size=20, face="bold"),
		 	legend.title=element_blank(),
		 	legend.key.size=unit(20, "points"))
invisible(dev.off())

png(paste(configOptions$prefix, "category_distribution_ts_byscore.png", sep="_"), width=3500, height=1500, res=300)
ggplot(subset(cancer.score.df, Score=="1" & Score.category=="Tumor Suppressor"), aes(x=Cancer, fill=Category)) + 
geom_histogram(binwidth=0.3, position="dodge") + 
facet_wrap( ~ Category, ncol=5) +
coord_flip() +
guides(fill=FALSE) + 
scale_fill_manual(values=score.colors) +
scale_y_continuous(breaks=c(2000, 6000, 10000), labels=c(2000, 6000, 10000)) +
labs(y="Number of genes", x="") + 
theme(axis.text.x=element_text(color="gray30", size=16, face="bold", vjust=0),
			axis.text.y=element_text(color="gray30", size=20, face="bold"),
		 	title=element_text(color="gray30", size=20, face="bold"),
		 	strip.text=element_text(color="gray30", size=20, face="bold"),
		 	legend.title=element_blank())
dev.off()

