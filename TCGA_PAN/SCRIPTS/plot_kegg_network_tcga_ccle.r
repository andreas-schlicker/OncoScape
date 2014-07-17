# Load packages
library(igraph)
library(plotrix)
library(ggplot2)

# Import scripts
source("/srv/nfs4/medoid-home/NKI/a.schlicker/R_SCRIPTS/biomart.r")

source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/utils.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/plotting.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/scoring.r")

source("/srv/nfs4/medoid-home/NKI/a.schlicker/randomwalk/networkplot.r")

# Load the data
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata")
results.tcga = results
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CCLE/20140416/prioritize_tcga_pancancer_allgenes_step2.rdata")
results.ccle = results

load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/NETWORK/KEGG/20140225/kegg_mat.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/NETWORK/KEGG/20140225/kegg_pathways.rdata")

# Build the KEGG adjacency matrix
#kegg.adj = matrix(0, ncol=ncol(kegg.mat), nrow=ncol(kegg.mat))
#colnames(kegg.adj) = colnames(kegg.mat)
#rownames(kegg.adj) = colnames(kegg.mat)
#for (i in 1:nrow(kegg.mat)) {
#	paths = names(which(kegg.mat[i, ] == 1))
#	kegg.adj[paths, paths] = 1
#}
#diag(kegg.adj) = 0


# Translate Entrez IDs to gene symbols
geneIds = translateIds(rownames(kegg.mat), filter="entrezgene")

# Filter the list of pathways
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

# Build the KEGG adjacency matrix
kegg.adj = matrix(0, ncol=length(kegg.list), nrow=length(kegg.list))
rownames(kegg.adj) = names(kegg.list)
colnames(kegg.adj) = names(kegg.list)
for (i in 1:length(kegg.list)) {
	for (j in (i+1):length(kegg.list)) {
		if (j > length(kegg.list)) {
			break
		}
		kegg.adj[i, j] = length(intersect(kegg.list[[i]], kegg.list[[j]]))
		kegg.adj[j, i] = length(intersect(kegg.list[[i]], kegg.list[[j]]))
	}
}
kegg.adj[which(kegg.adj < 3)] = 0
kegg.adj[which(kegg.adj > 0)] = 1


# Save the gene-based version of the results
results.tcga.gene = results.tcga
results.ccle.gene = results.ccle
# Create the pathway-based version of the results
results.tcga = lapply(results.tcga, function(x) { list(prioritize.combined=scoreGenesets(x$prioritize.combined, kegg.list)) })
results.ccle = lapply(results.ccle, function(x) { list(prioritize.combined=scoreGenesets(x$prioritize.combined, kegg.list)) })

# For the scatter plot
plotting.df = data.frame()
label.df = data.frame()

# Create a KEGG network plot for each cancer type
for (cancType in names(results.ccle)) {
	#dir.create(cancType)
	
	# Get the scores and thresholds
	tempScores = results.tcga[[cancType]]$prioritize.combined[, "combined.score"]
	thresholds = quantile(tempScores, probs=c(0.1, 0.9))
	# Only keep pathways with scores that are more extreme than the thresholds
	retain = unique(names(tempScores)[which(tempScores < thresholds[1] | tempScores > thresholds[2])])
	
	# Get the scores and thresholds
	tempScores.ccle = results.ccle[[cancType]]$prioritize.combined[, "combined.score"]
	thresholds.ccle = quantile(tempScores.ccle, probs=c(0.1, 0.9))
	# Only keep pathways with scores that are more extreme than the thresholds
	retain = union(retain, unique(names(tempScores.ccle)[which(tempScores.ccle < thresholds.ccle[1] | tempScores.ccle > thresholds.ccle[2])]))
	
	tmpAdj = kegg.adj[retain, retain]
	rownames(tmpAdj) = gsub(" signaling$", " sig.", gsub(" pathway", "", kegg.pathways[rownames(tmpAdj), 2])) 
	colnames(tmpAdj) = gsub(" signaling$", " sig.", gsub(" pathway", "", kegg.pathways[colnames(tmpAdj), 2]))
	
	# Create the graph
	g = createGraph(tmpAdj, tempScores[retain], mode="undirected", weighted=FALSE)
	# Plot sample graphs for each expression situation
	g = nodeColsByWeight(g, low="#034b87", mid="gray80", high="#880000")
	V(g)$size = 30
	g = setLayout(g, layout.fruchterman.reingold(g, params=list(area=vcount(g)^6)))
	png(paste(cancType, "_tcga_kegg_combined.score.png", sep=""), width=2000, height=2000, res=300)
	#tkplot(g, edge.width=5, vertex.label.color="gray20", vertex.label.family="arial",
	#		 vertex.label.font=2, vertex.label.dist=1, vertex.label.cex=1.1)
	plot(g, edge.width=5, vertex.label.color="gray20", vertex.label.family="arial",
			 vertex.label.font=2, vertex.label.dist=1, vertex.label.cex=1.1)
	color.legend(-1.5, 1.3, -1, 1.4, 
							 legend=signif(range(V(g)$weight), digits=2), 
							 rect.col=unique(V(g)$color[order(V(g)$weight)]), 
							 gradient="x", align="rb", cex=1.3, font=2)
	dev.off()
	
	V(g)$weight = tempScores.ccle[retain]
	g = nodeColsByWeight(g, low="#034b87", mid="gray80", high="#880000")
	png(paste(cancType, "_ccle_kegg_combined.score.png", sep=""), width=2000, height=2000, res=300)
	#tkplot(g, edge.width=5, vertex.label.color="gray20", vertex.label.family="arial",
	#		 vertex.label.font=2, vertex.label.dist=1, vertex.label.cex=1.1)
	plot(g, edge.width=5, vertex.label.color="gray20", vertex.label.family="arial",
			 vertex.label.font=2, vertex.label.dist=1, vertex.label.cex=1.1)
	color.legend(-1.5, 1.3, -1, 1.4, 
							 legend=signif(range(V(g)$weight), digits=2), 
							 rect.col=unique(V(g)$color[order(V(g)$weight)]), 
							 gradient="x", align="rb", cex=1.3, font=2)
	dev.off()

	generatePathview(results.tcga.gene, retain, cancers=cancType, scores="combined.score", 
									 out.dir=cancType, out.suffix="combined.score",
									 kegg.dir="/srv/nfs4/medoid-bulk/NKI/a.schlicker/EXTERNAL_DB/KEGG", multi.state=FALSE,
									 colors=list(low="#034b87", mid="gray98", high="#880000"))
									 
	generatePathview(results.ccle.gene, retain, cancers=cancType, scores="combined.score", 
									 out.dir=cancType, out.suffix="ccle_combined.score",
									 kegg.dir="/srv/nfs4/medoid-bulk/NKI/a.schlicker/EXTERNAL_DB/KEGG", multi.state=FALSE,
									 colors=list(low="#034b87", mid="gray98", high="#880000"))
									 
	plotting.df = rbind(plotting.df,
											data.frame(CCLE=tempScores.ccle,
													 			 TCGA=tempScores[names(tempScores.ccle)],
													 			 Pathway=names(tempScores.ccle),
													 			 Cancer=cancType))
	label.df = rbind(label.df,
									 data.frame(Cancer=cancType,
									 						label=paste("cor=", signif(cor(tempScores.ccle, tempScores[names(tempScores.ccle)], method="spearman"), digits=2), sep=""),
									 						x=-0.4, y=0.4))
}
plotting.df$Cancer = factor(plotting.df$Cancer, levels=sort(unique(as.character(plotting.df$Cancer))))
label.df$Cancer = factor(label.df$Cancer, levels=sort(unique(as.character(label.df$Cancer))))

png("kegg_tcga_ccle_scatter.png", width=5000, height=3000, res=300)
ggplot(plotting.df, aes(x=TCGA, y=CCLE)) +
geom_point(size=4) +
geom_text(data=label.df, aes(x=x, y=y, label=label), hjust=0, fontface="bold", size=8, show_guide=FALSE) +
labs(x="Avg. combined score TCGA", y="Avg. combined score CCLE") + 
facet_wrap(~ Cancer, ncol=5) +
theme(axis.text=element_text(size=20, face="bold", color="gray20"),
			axis.title.x=element_text(size=20, face="bold", color="gray20", vjust=0),
			axis.title.y=element_text(size=20, face="bold", color="gray20", vjust=0.3),
			strip.text=element_text(size=20, face="bold", color="gray20"))
dev.off()

