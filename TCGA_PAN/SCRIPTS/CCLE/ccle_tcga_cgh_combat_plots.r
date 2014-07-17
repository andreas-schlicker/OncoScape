library(fastcluster)
library(NMF)
library(gplots)

rownames(sample.ann) = c(mappedCcle, 
												 unlist(lapply(names(tcga.cna.combat), function(x) { colnames(tcga.cna.combat[[x]][[1]]) })), 
												 unlist(lapply(names(tcga.mn.cna.combat), function(x) { colnames(tcga.mn.cna.combat[[x]][[1]]) })))

sample.ann = subset(sample.ann, tissue != "LAML")

for (cancType in unique(sample.ann[, 3])) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]
	
	temp = cbind(ccle.cna.combat[[cancType]][[1]],
							 tcga.cna.combat[[cancType]][[1]],
							 tcga.mn.cna.combat[[cancType]][[1]])
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ])), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]))), method="complete"))
	
	ext = max(abs(min(temp[vars, ])), max(temp[vars, ]))
	
	png(paste("tcga_ccle_cna_", cancType, "_combat.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", mid="white", high="yellow"),
         	 breaks=seq(-1*ext, ext, length.out=50),
         	 Colv=colv,
         	 Rowv=rowv,
         	 distfun="pearson",
         	 annCol=sample.ann[colnames(temp), ],
         	 labRow=NA,
         	 labCol=NA)
  dev.off()
}

for (cancType in unique(sample.ann[, 3])) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]
	
	temp = cbind(ccle.cna.combat[[cancType]][[1]],
							 tcga.mn.cna.combat[[cancType]][[1]])
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ])), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]))), method="complete"))
	
	ext = max(abs(min(temp[vars, ])), max(temp[vars, ]))
	
	png(paste("tcga_ccle_cna_", cancType, "_combat_nots.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", mid="white", high="yellow"),
         	 breaks=seq(-1*ext, ext, length.out=50),
         	 Colv=colv,
         	 Rowv=rowv,
         	 distfun="pearson",
         	 annCol=sample.ann[colnames(temp), ],
         	 labRow=NA,
         	 labCol=NA)
  dev.off()
}

for (cancType in setdiff(unique(sample.ann[, 3]), "READ"))) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]
	
	commonGenes = intersect(rownames(ccle.cna[[cancType]][[1]]),
													intersect(rownames(tcga.cna[[cancType]][[1]]),
																		rownames(tcga.mn.cna[[cancType]][[1]])))
	
	temp = cbind(impute.knn(ccle.cna[[cancType]][[1]][commonGenes, ])$data, 
							 impute.knn(tcga.cna[[cancType]][[1]][commonGenes, ])$data, 
							 impute.knn(tcga.mn.cna[[cancType]][[1]][commonGenes, ])$data)
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ], use="pairwise.complete.obs")), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]), use="pairwise.complete.obs")), method="complete"))
	
	ext = max(abs(min(temp[vars, ])), max(temp[vars, ]))
	
	png(paste("tcga_ccle_cna_", cancType, "_nocombat.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", mid="white", high="yellow"),
         	 breaks=seq(-1*ext, ext, length.out=50),
         	 Colv=colv,
         	 Rowv=rowv,
         	 distfun="pearson",
         	 annCol=sample.ann[colnames(temp), ],
         	 labRow=NA,
         	 labCol=NA)
  dev.off()
}

for (cancType in setdiff(unique(sample.ann[, 3]), "READ")) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]
	
	commonGenes = intersect(rownames(ccle.cna[[cancType]][[1]]),
													intersect(rownames(tcga.cna[[cancType]][[1]]),
																		rownames(tcga.mn.cna[[cancType]][[1]])))
	
	temp = cbind(impute.knn(ccle.cna[[cancType]][[1]][commonGenes, ])$data, 
							 impute.knn(tcga.mn.cna[[cancType]][[1]][commonGenes, ])$data)
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ], use="pairwise.complete.obs")), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]), use="pairwise.complete.obs")), method="complete"))
	
	ext = max(abs(min(temp[vars, ])), max(temp[vars, ]))
	
	png(paste("tcga_ccle_cna_", cancType, "_nocombat_nots.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", mid="white", high="yellow"),
         	 breaks=seq(-1*ext, ext, length.out=50),
         	 Colv=colv,
         	 Rowv=rowv,
         	 distfun="pearson",
         	 annCol=sample.ann[colnames(temp), ],
         	 labRow=NA,
         	 labCol=NA)
  dev.off()
}

