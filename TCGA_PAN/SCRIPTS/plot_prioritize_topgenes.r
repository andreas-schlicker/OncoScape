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


for (n in c("og", "ts")) {
	result.df = heatmapDataframe(results, scores=list(combined=paste(n, ".score", sep=""),
																										Meth=paste(n, ".methylation", sep=""), 
																										CNA=paste(n, ".cna", sep=""), 
																										Mut=paste(n, ".mutations", sep=""), 
																										shRNA=paste(n, ".achilles", sep=""), 
																										Expr=paste(n, ".exprs", sep="")))
	
	gene.order = names(sort(unlist(lapply(split(subset(result.df, score.type=="combined")[, "score"], 
																							subset(result.df, score.type=="combined")[, "gene"]), 
																				sum, 
																				na.rm=TRUE))))
	topgenes = gene.order[max(1, (length(gene.order)-configOptions$notopgenes-1)):length(gene.order)]
	result.df$gene = factor(result.df$gene, levels=gene.order)
	result.df$score.type = factor(result.df$score.type, levels=c("CNA", "Expr", "Meth", "Mut", "shRNA", "combined"))
	result.df$cancer = factor(result.df$cancer, levels=sort(unique(as.character(result.df$cancer))))
	
	result.df[, 2] = as.character(result.df[, 2])
	result.df[which(!is.na(result.df[, 2]) & result.df[, 2] == "1"), 2] = as.character(result.df[which(!is.na(result.df[, 2]) & result.df[, 2] == "1"), 3])
	result.df[which(is.na(result.df[, 2]) | result.df[, 2] == "0"), 2] = "NONE"
	
	png(paste(configOptions$prefix, "_top_", n, "_genes_details.png", sep=""), width=4000, height=2500, res=300)
	print(
	ggplot(subset(result.df, score.type != "combined" & gene %in% topgenes), aes(x=score.type, y=gene)) + 
	geom_tile(aes(fill=score), color="white", size=0.7) +
	scale_fill_manual(values=c(NONE="white", CNA="#888888", Expr="#E69F00", Meth="#56B4E9", Mut="#009E73", shRNA="#F0E442"), 
										breaks=c("CNA", "Expr", "Meth", "Mut", "shRNA")) +
	labs(x="", y="") +
	facet_grid(.~cancer) + 
	theme(panel.background=element_rect(color="white", fill="white"),
				panel.margin=unit(10, "points"),
				axis.ticks=element_blank(), 
				axis.text.x=element_blank(),
				axis.text.y=element_text(color="gray30", size=16, face="bold"),
				axis.title.x=element_text(color="gray30", size=16, face="bold"),
				strip.text.x=element_text(color="gray30", size=16, face="bold"),
				legend.text=element_text(color="gray30", size=16, face="bold"),
				legend.title=element_blank(),
				legend.position="bottom")
	)
	dev.off()
	
	summaryTable = table(subset(result.df, score.type != "combined" & gene %in% topgenes)[, c("gene", "score")])[topgenes, ]
	summaryDf = data.frame()
	for (g in rownames(summaryTable)) {
		summaryDf = rbind(summaryDf, 
											data.frame(gene=g,
																 score=c("CNA", "Expr", "Meth", "Mut", "shRNA"),
																 value=summaryTable[g, c("CNA", "Expr", "Meth", "Mut", "shRNA")]))
	}

	summaryDf$gene = factor(summaryDf$gene, levels=rev(topgenes))

	png(paste(configOptions$prefix, "_top_", n, "_genes_details_marginalhist.png", sep=""), width=500, height=2500, res=300)
	print(
	ggplot(summaryDf, aes(x=score, y=value, fill=score)) + 
	geom_histogram(stat="identity") + 
	coord_flip() +
	guides(fill=FALSE) +
	scale_fill_manual(values=c(NONE="white", CNA="#888888", Expr="#E69F00", Meth="#56B4E9", Mut="#009E73", shRNA="#F0E442"), 
											breaks=c("CNA", "Expr", "Meth", "Mut", "shRNA")) +
	scale_y_continuous(breaks=seq(2, 8, by=2)) +
	facet_grid(gene ~ .) + 
	theme(panel.background=element_rect(color="white", fill="white"),
				panel.grid=element_blank(),
				strip.text=element_blank(),
				axis.title=element_blank(),
				axis.text.x=element_text(size=16, face="bold", color="gray30"),
				axis.text.y=element_blank(),
				axis.ticks.y=element_blank())
	)
	dev.off()
}


