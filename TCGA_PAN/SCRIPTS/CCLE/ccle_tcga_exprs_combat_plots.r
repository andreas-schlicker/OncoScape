library(fastcluster)
library(NMF)
library(gplots)

mappedCcle = intersect(as.character(cls.df[, "CL"]), colnames(cl.ccle.exprs))
												 
sample.ann = rbind(data.frame(batch="ccle", cancer="cancer", tissue=cls.df[mappedCcle, 1]),
									 data.frame(batch="tcga", 
										 					cancer="cancer", 
										 					tissue=unlist(lapply(names(tcga.exprs), function(x) { rep(x, times=ncol(tcga.exprs[[x]][[1]])) }))),
									 data.frame(batch="tcga", 
										 					cancer="normal", 
										 					tissue=unlist(lapply(names(tcga.mn.exprs), function(x) { rep(x, times=ncol(tcga.mn.exprs[[x]][[1]])) }))))

rownames(sample.ann) = c(mappedCcle, 
												 unlist(lapply(names(tcga.exprs), function(x) { colnames(tcga.exprs[[x]][[1]]) })), 
												 unlist(lapply(names(tcga.mn.exprs), function(x) { colnames(tcga.mn.exprs[[x]][[1]]) })))

sample.ann = subset(sample.ann, tissue != "LAML")

for (cancType in unique(sample.ann[, 3])) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]
	
	temp = cbind(ccle.exprs.combat[[cancType]][[1]],
							 tcga.exprs.combat[[cancType]][[1]],
							 tcga.mn.exprs.combat[[cancType]][[1]])
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ])), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]))), method="complete"))
	
	png(paste("tcga_ccle_exprs_", cancType, "_combat.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", high="yellow"),
         	 breaks=seq(min(temp[vars, ]), 
         	 						max(temp[vars, ]), 
         	 						length.out=50),
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
	
	temp = cbind(ccle.exprs.combat[[cancType]][[1]],
							 tcga.mn.exprs.combat[[cancType]][[1]])
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ])), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]))), method="complete"))
	
	png(paste("tcga_ccle_exprs_", cancType, "_combat_nots.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", high="yellow"),
         	 breaks=seq(min(temp[vars, ]), 
         	 						max(temp[vars, ]), 
         	 						length.out=50),
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
	
	commonGenes = intersect(rownames(ccle.exprs[[cancType]][[1]]),
													intersect(rownames(tcga.exprs[[cancType]][[1]]),
																		rownames(tcga.mn.exprs[[cancType]][[1]])))
	
	temp = cbind(impute.knn(ccle.exprs[[cancType]][[1]][commonGenes, ])$data, 
							 impute.knn(tcga.exprs[[cancType]][[1]][commonGenes, ])$data, 
							 impute.knn(tcga.mn.exprs[[cancType]][[1]][commonGenes, ])$data)
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ], use="pairwise.complete.obs")), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]), use="pairwise.complete.obs")), method="complete"))
	
	png(paste("tcga_ccle_exprs_", cancType, "_nocombat.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", high="yellow"),
         	 breaks=seq(min(temp[vars, ]), 
         	 						max(temp[vars, ]), 
         	 						length.out=50),
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
	
	commonGenes = intersect(rownames(ccle.exprs[[cancType]][[1]]),
													rownames(tcga.mn.exprs[[cancType]][[1]]))
	
	temp = cbind(impute.knn(ccle.cna[[cancType]][[1]][commonGenes, ])$data, 
							 impute.knn(tcga.mn.cna[[cancType]][[1]][commonGenes, ])$data)
	
	vars = names(sort(apply(temp, 1, IQR, na.rm=TRUE), decreasing=TRUE))[1:5000]

	colv=as.dendrogram(hclust(as.dist(1-cor(temp[vars, ], use="pairwise.complete.obs")), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[vars, ]), use="pairwise.complete.obs")), method="complete"))
	
	png(paste("tcga_ccle_exprs_", cancType, "_nocombat_nots.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(temp[vars, ],
         	 scale="none",
         	 color=colorpanel(49, low="blue", high="yellow"),
         	 breaks=seq(min(temp[vars, ]), 
         	 						max(temp[vars, ]), 
         	 						length.out=50),
         	 Colv=colv,
         	 Rowv=rowv,
         	 distfun="pearson",
         	 annCol=sample.ann[colnames(temp), ],
         	 labRow=NA,
         	 labCol=NA)
  dev.off()
}

