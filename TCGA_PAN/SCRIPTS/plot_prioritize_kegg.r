#!/usr/bin/Rscript

library(stringr)
library(ggplot2)
library(grid)

source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/utils.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/plotting.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/scoring.r")

source("/srv/nfs4/medoid-home/NKI/a.schlicker/R_SCRIPTS/biomart.r")

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


# Read the KEGG mapping and prepare the matrices
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/NETWORK/KEGG/20140225/kegg_mat.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/NETWORK/KEGG/20140225/kegg_pathways.rdata")

# Map Entrez IDs to gene symbols
geneIds = translateIds(rownames(kegg.mat), filter="entrezgene")

# Take out only those pathways starting with hsa04 and remove genes that don't map there
pathways = c("hsa04514", "hsa04010", "hsa04012", "hsa04014", "hsa04015", "hsa04062", "hsa04066", "hsa04150", "hsa04151", 
						 "hsa04210", "hsa04370", "hsa04510", "hsa04530", "hsa04620", "hsa04630", "hsa04668", "hsa04910", "hsa04310", 
						 "hsa04120", "hsa04810", "hsa04020", "hsa04520", "hsa04110", "hsa04115", "hsa04390", "hsa04911", "hsa04350", 
						 "hsa04512", "hsa04060", "hsa04621", "hsa04070", "hsa04064", "hsa04330", "hsa04340")
kegg.mat = data.matrix(kegg.mat[, intersect(colnames(kegg.mat), pathways)])
kegg.mat = kegg.mat[apply(kegg.mat, 1, function(x) { sum(x) > 0 }), ]
# Replace Entrez ID with gene symbol
for (i in 1:nrow(kegg.mat)) {
  rownames(kegg.mat)[i] = geneIds[which(geneIds[, "entrezgene"] == as.character(rownames(kegg.mat)[i])), "hgnc_symbol"][1]
}
# And make a list out of the pathways
kegg.list = lapply(colnames(kegg.mat), function(x) { rownames(kegg.mat)[which(kegg.mat[, x] == 1)] })
names(kegg.list) = colnames(kegg.mat)

results.gene = results

results = lapply(results, function(x) { list(prioritize.combined=scoreGenesets(x$prioritize.combined, kegg.list)) })

MAXGENES = min(1000, nrow(results[[1]]$prioritize.combined))

result.df = heatmapDataframe(results)
result.df[which(result.df[, "score.type"] == "TS"), "score"] = -1 * result.df[which(result.df[, "score.type"] == "TS"), "score"]
result.df[, "gene"] = kegg.pathways[as.character(result.df[, "gene"]), 2]



## Need to add indexes of topgenes as parameters
parameters = list()

gene.order = getOrder(result.df, "TS", "gene")
topgenes = gene.order[1:min(MAXGENES, length(gene.order))]
cancer.order = getOrder(result.df, "TS", "cancer")
parameters[["ts1"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=18)), score.type="TS", 
													 color.low="#034b87", color.high="gray98",
													 ylab="", title="Tumor suppressor score")
parameters[["ts2"]] = list(data=result.df, prefix=configOptions$prefix, filename="ts_affected_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=18)), score.type="TS.affected", 
													 color.low="gray98", color.high="#034b87",
													 ylab="", title="% affected samples TS")
													 
gene.order = getOrder(result.df, "OG", "gene")
topgenes = gene.order[max(1, (length(gene.order)-MAXGENES-1)):length(gene.order)]
cancer.order = getOrder(result.df, "OG", "cancer")
parameters[["og1"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=18)), score.type="OG",
													 color.low="gray98", color.high="#880000", 
													 ylab="", title="Oncogene score")
parameters[["og2"]] = list(data=result.df, prefix=configOptions$prefix, filename="og_affected_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=18)), score.type="OG.affected",
													 color.low="gray98", color.high="#880000", 
													 ylab="", title="% affected samples OG")
													 
gene.order = getOrder(result.df, "Combined", "gene")
topgenes= if (MAXGENES == length(gene.order)) {	
						topgenes = gene.order 
					} else { 
						topgenes = c(gene.order[1:ceiling(MAXGENES/2)], gene.order[(length(gene.order)-(floor(MAXGENES/2)-1)):length(gene.order)]) 
					}
cancer.order = getOrder(result.df, "Combined", "cancer")
parameters[["cs1"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=18)), score.type="Combined",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Combined score")
parameters[["cs2"]] = list(data=result.df, prefix=configOptions$prefix, filename="combined_affected_heatmap.png", gene.order=gene.order,
													 topgenes=topgenes, cancer.order=cancer.order,
													 yaxis=theme(axis.text.y=element_text(color="gray30", face="bold", size=18)), score.type="Combined.affected",
													 color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="% affected samples Combined score")

topgenes= c(gene.order[1:configOptions$notopgenes], gene.order[(length(gene.order)-configOptions$notopgenes+1):length(gene.order)])

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

png(paste(configOptions$prefix, "score_histogram_by_pathway.png", sep="_"), width=4000, height=3000, res=300)
print(getDistPlot(subset(result.df, score.type %in% c("OG", "TS", "Combined") & gene %in% topgenes), facets="~ gene",
									plot.type="histogram", ncol=4, title="Score distribution by pathway", x="Score", y="Number of cancers",
									strip.theme=theme(strip.text=element_text(size=18))))
invisible(dev.off())

png(paste(configOptions$prefix, "score_density_by_pathway.png", sep="_"), width=4000, height=3000, res=300)
print(getDistPlot(subset(result.df, score.type %in% c("OG", "TS", "Combined") & gene %in% topgenes), facets="~ gene",
									plot.type="density", ncol=4, title="Score distribution by pathway", x="Score", y="Density",
									strip.theme=theme(strip.text=element_text(size=18))))
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
																		 		 Score.category=toupper(st),
																		 		 Score=results[[cancType]]$prioritize.combined[, paste(st, score, sep=".")]))
		}
	}
}

png(paste(configOptions$prefix, "category_distribution.png", sep="_"), width=3000, height=3000, res=300)
ggplot(cancer.score.df, aes(x=Score, fill=Category)) + 
geom_histogram(binwidth=0.1, position="dodge") + 
facet_grid(Cancer ~ Score.category) +
scale_fill_manual(values=score.colors2) + 
labs(title="Score distribution by gene", y="Number of genes", x="Score") + 
theme(axis.ticks=element_blank(), 
		 	axis.text.x=element_text(colour="gray30", face="bold"),
		 	axis.text.y=element_text(colour="gray30", face="bold"))
invisible(dev.off)


## Plots of pathway maps
out.dir = dirname(configOptions$prefix)
paths = c(topgenes[1:5], topgenes[(length(topgenes)-4):length(topgenes)])
paths = kegg.pathways[which(kegg.pathways[, 2] %in% paths), 1]
for (score in c("combined.score", "og.score", "ts.score")) {
	multiplier = 1
	if (score == "ts.score") {
		multiplier = -1
	}
	
	generatePathview(results.gene, paths, scores=score, out.dir=out.dir, out.suffix=score,
									 kegg.dir="/srv/nfs4/medoid-bulk/NKI/a.schlicker/EXTERNAL_DB/KEGG", multi.state=TRUE, multiplier=multiplier,
									 colors=list(low="#034b87", mid="gray98", high="#880000"))
}

