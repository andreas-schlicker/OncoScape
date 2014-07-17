#!/usr/bin/Rscript

setwd("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/KEGG_PATHWAY_HEATMAPS")

outputHeatmap = function(params) {
	p = ggplot(subset(params$data, score.type==params$score.type), aes(x=set, y=gene)) + 
			 			 geom_tile(aes(fill=score), color = "white") +
			 			 facet_grid(. ~ cancer) +
			 			 labs(title=params$title, x="", y="") +
			 			 theme(panel.background=element_rect(color="white", fill="white"),
				  	 			 axis.ticks=element_blank(), 
				  	 			 axis.text.x=element_text(color="gray30", face="bold", size=20, angle=90),
				  	 			 axis.text.y=element_text(color="gray30", face="bold", size=18),
				  	 			 strip.text=element_text(color="gray30", face="bold", size=20),
				  	 			 legend.text=element_text(color="gray30", face="bold", size=20),
				  	 			 legend.title=element_blank(),
				  	 			 legend.position="top",
				  	 			 legend.key.width=unit(1.5, "cm"))
	if (!is.null(params$color.mid)) {
		p = p + scale_fill_gradient2(low=params$color.low, mid=params$color.mid, high=params$color.high)
	} else {
		p = p + scale_fill_gradient(low=params$color.low, high=params$color.high)
	}
	
	png(params$filename, width=5000, height=3000, res=300)
	print(p)
	invisible(dev.off())
}

library(stringr)
library(ggplot2)
library(grid)

source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/utils.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/plotting.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/scoring.r")

source("/srv/nfs4/medoid-home/NKI/a.schlicker/R_SCRIPTS/biomart.r")


# Load prioritization results
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata")
results.tcga = results
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CCLE/20140416/prioritize_tcga_pancancer_allgenes_step2.rdata")
results.ccle = results

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


results.tcga.kegg = lapply(results.tcga, function(x) { list(prioritize.combined=scoreGenesets(x$prioritize.combined, kegg.list)) })
results.ccle.kegg = lapply(results.ccle, function(x) { list(prioritize.combined=scoreGenesets(x$prioritize.combined, kegg.list)) })

result.df.tcga = cbind(heatmapDataframe(results.tcga.kegg), set="TCGA")
#result.df.tcga = paste(as.character(result.df.tcga$cancer), "tcga", sep="_")
result.df.ccle = cbind(heatmapDataframe(results.ccle.kegg), set="CCLE")
#result.df.ccle$cancer = paste(as.character(result.df.ccle$cancer), "ccle", sep="_")

result.df = rbind(result.df.tcga, result.df.ccle)
result.df$cancer = factor(result.df$cancer, levels=sort(levels(result.df$cancer)))
result.df[which(result.df[, "score.type"] == "TS"), "score"] = -1 * result.df[which(result.df[, "score.type"] == "TS"), "score"]
result.df$gene = gsub(" pathway", "", kegg.pathways[as.character(result.df$gene), 2])
result.df$gene = factor(result.df$gene, levels=sort(unique(result.df$gene), decreasing=TRUE))


parameters = list()

gene.order = getOrder(result.df, "TS", "gene")
parameters[["ts1"]] = list(data=result.df, filename="kegg_tcga_ccle_ts_heatmap.png",
													 score.type="TS", color.low="#034b87", color.high="gray98",
													 ylab="", title="Tumor suppressor score",
													 xaxis=theme(axis.text.x=element_text(color="gray30", face="bold", size=18, angle=90)))													 
													 
gene.order = getOrder(result.df, "OG", "gene")
parameters[["og1"]] = list(data=result.df, filename="kegg_tcga_ccle_og_heatmap.png",
													 score.type="OG", color.low="gray98", color.high="#880000", 
													 ylab="", title="Oncogene score",
													 xaxis=theme(axis.text.x=element_text(color="gray30", face="bold", size=18, angle=90)))
													 
													 
gene.order = getOrder(result.df, "Combined", "gene")
#cancer.order = getOrder(result.df, "Combined", "cancer")
parameters[["cs1"]] = list(data=result.df, filename="kegg_tcga_ccle_combined_heatmap.png",
													 score.type="Combined", color.low="#034b87", color.mid="gray98", color.high="#880000",
													 ylab="", title="Combined score",
													 xaxis=theme(axis.text.x=element_text(color="gray30", face="bold", size=18, angle=90)))
													 				 
invisible(lapply(parameters, outputHeatmap))


