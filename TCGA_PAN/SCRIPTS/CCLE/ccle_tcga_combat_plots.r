library(fastcluster)
library(NMF)
library(gplots)

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

complete.exprs = matrix(NA, nrow=length(commonGenes), ncol=length(mappedCcle) +
																														 sum(unlist(lapply(tcga.exprs, function(x) { ncol(x[[1]]) }))) + 
																														 sum(unlist(lapply(tcga.mn.exprs, function(x) { ncol(x[[1]]) }))))
	rownames(complete.exprs) = commonGenes
	colnames(complete.exprs) = rep("", times=ncol(complete.exprs))

	i = 1
	last = i + length(mappedCcle) - 1
	colnames(complete.exprs)[i:last] = mappedCcle
	complete.exprs[, i:last] = cl.ccle.exprs[commonGenes, mappedCcle] - apply(cl.ccle.exprs[commonGenes, mappedCcle], 1, mean)

	for (j in names(tcga.exprs)) {
		i = last + 1
		last = i + ncol(tcga.exprs[[j]][[1]]) - 1
		colnames(complete.exprs)[i:last] = colnames(tcga.exprs[[j]][[1]])

		temp = tcga.exprs[[j]][[1]][commonGenes, ]
		if (length(which(is.na(temp))) > 0) {
			temp = impute.knn(temp)$data
		}
	
		if (min(temp, na.rm=TRUE) >= 0) {	
			temp = log2(temp + 1) - apply(log2(temp + 1), 1, mean)
		}

		complete.exprs[, i:last] = temp
	}

	for (j in names(tcga.mn.exprs)) {
		i = last + 1
		last = i + ncol(tcga.mn.exprs[[j]][[1]]) - 1
		colnames(complete.exprs)[i:last] = colnames(tcga.mn.exprs[[j]][[1]])
		
		temp = tcga.mn.exprs[[j]][[1]][commonGenes, ]
		if (length(which(is.na(temp))) > 0) {
			temp = impute.knn(temp)$data
		}
		
		if (min(temp, na.rm=TRUE) >= 0) {	
			temp = log2(temp + 1) - apply(log2(temp + 1), 1, mean)
		}
	
		complete.exprs[, i:last] = temp	
	}

												 
vars = apply(complete.exprs, 1, IQR, na.rm=TRUE)
												 
colv=as.dendrogram(hclust(as.dist(1-cor(complete.exprs[names(sort(vars, decreasing=TRUE)[1:1000]),])), method="complete"))
rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.exprs[names(sort(vars, decreasing=TRUE)[1:1000]),]))), method="complete"))
         
png("tcga_ccla_exprs_new.png", width=4000, height=3000, res=300)
aheatmap(complete.exprs[names(sort(vars, decreasing=TRUE)[1:1000]),],
         scale="none",
         color=colorpanel(49, low="blue", mid="white", high="yellow"),
         breaks=seq(-8, 5, length.out=50),
         Colv=colv,
         Rowv=rowv,
         distfun="pearson",
         annColors=list(tissue=c("#000000", "#555555", "#999999", "#CCCCCC", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")),
         annCol=sample.ann,
         labRow=NA,
         labCol=NA)
dev.off()

colv=as.dendrogram(hclust(as.dist(1-cor(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:1000]),])), method="complete"))
rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:1000]),]))), method="complete"))
         
png("tcga_ccla_exprs_combat.png", width=4000, height=3000, res=300)
aheatmap(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:1000]),],
         scale="none",
         color=colorpanel(49, low="blue", mid="white", high="yellow"),
         breaks=seq(-8, 5, length.out=50),
         Colv=colv,
         Rowv=rowv,
         distfun="pearson",
         annCol=sample.ann,
         labRow=NA,
         labCol=NA)
dev.off()

												 
												 
for (cancType in unique(sample.ann[, 3])) {
	#samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType & 
																			 (sample.ann[, "batch"] == "ccle" | (sample.ann[, "batch"] == "tcga" & sample.ann[, "cancer"] == "normal")) )]
	
	vars = apply(complete.exprs[, samples], 1, IQR, na.rm=TRUE)

	colv=as.dendrogram(hclust(as.dist(1-cor(complete.exprs[names(sort(vars, decreasing=TRUE)[1:5000]), samples])), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.exprs[names(sort(vars, decreasing=TRUE)[1:5000]), samples]))), method="complete"))
	
	png(paste("tcga_ccla_exprs_", cancType, ".png", sep=""), width=4000, height=3000, res=300)
	aheatmap(complete.exprs[names(sort(vars, decreasing=TRUE)[1:5000]), samples],
         	 scale="none",
         	 color=colorpanel(49, low="blue", high="yellow"),
         	 breaks=seq(min(complete.exprs[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
         	 						max(complete.exprs[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
         	 						length.out=50),
         	 Colv=colv,
         	 Rowv=rowv,
         	 distfun="pearson",
         	 annCol=sample.ann[samples, ],
         	 labRow=NA,
         	 labCol=NA)
  dev.off()


  vars = apply(complete.exprs.combat[, samples], 1, IQR, na.rm=TRUE)

  colv=as.dendrogram(hclust(as.dist(1-cor(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples])), method="complete"))
  rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples]))), method="complete"))
         
  png(paste("tcga_ccla_exprs_", cancType, "_combat.png", sep=""), width=4000, height=3000, res=300)
  aheatmap(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples],
  	       scale="none",
   	       color=colorpanel(49, low="blue", high="yellow"),
   	       breaks=seq(min(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						max(complete.exprs.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						length.out=50),
   	       Colv=colv,
   	       Rowv=rowv,
   	       distfun="pearson",
   	       annCol=sample.ann[samples, ],
   	       labRow=NA,
   	       labCol=NA)
  dev.off()
}



for (cancType in setdiff(unique(sample.ann[, 3]), c("LAML", "HNSC"))) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]

	temp = complete.exprs[, samples]
	mod = model.matrix(~as.factor(cancer), data=data.frame(sample.ann[samples, c("cancer", "tissue")]))								 
	# Normalize for TCGA vs CCLE
	temp = ComBat(temp, batch=sample.ann[, "batch"], mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

  vars = apply(temp, 1, IQR, na.rm=TRUE)

  colv=as.dendrogram(hclust(as.dist(1-cor(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples])), method="complete"))
  rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples]))), method="complete"))
         
  png(paste("tcga_ccla_exprs_", cancType, "_sepcombat.png", sep=""), width=4000, height=3000, res=300)
  aheatmap(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples],
  	       scale="none",
   	       color=colorpanel(49, low="blue", high="yellow"),
   	       breaks=seq(min(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						max(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						length.out=50),
   	       Colv=colv,
   	       Rowv=rowv,
   	       distfun="pearson",
   	       annCol=sample.ann[samples, ],
   	       labRow=NA,
   	       labCol=NA)
  dev.off()
}










rownames(sample.ann) = c(mappedCcle, 
												 unlist(lapply(names(tcga.cna), function(x) { colnames(tcga.cna[[x]][[1]]) })), 
												 unlist(lapply(names(tcga.mn.cna), function(x) { colnames(tcga.mn.cna[[x]][[1]]) })))

vars = apply(complete.cna, 1, IQR, na.rm=TRUE)
												 
colv=as.dendrogram(hclust(as.dist(1-cor(complete.cna[names(sort(vars, decreasing=TRUE)[1:1000]),])), method="complete"))
rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.cna[names(sort(vars, decreasing=TRUE)[1:1000]),]))), method="complete"))
         
png("tcga_ccla_cna.png", width=4000, height=3000, res=300)
aheatmap(complete.cna[names(sort(vars, decreasing=TRUE)[1:1000]),],
         scale="none",
         color=colorpanel(49, low="blue", mid="white", high="yellow"),
         breaks=seq(-8, 5, length.out=50),
         Colv=colv,
         Rowv=rowv,
         distfun="pearson",
         annCol=sample.ann,
         labRow=NA,
         labCol=NA)
dev.off()

colv=as.dendrogram(hclust(as.dist(1-cor(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:1000]),])), method="complete"))
rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:1000]),]))), method="complete"))
         
png("tcga_ccla_cna_combat.png", width=4000, height=3000, res=300)
aheatmap(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:1000]),],
         scale="none",
         color=colorpanel(49, low="blue", mid="white", high="yellow"),
         breaks=seq(-8, 5, length.out=50),
         Colv=colv,
         Rowv=rowv,
         distfun="pearson",
         annCol=sample.ann,
         labRow=NA,
         labCol=NA)
dev.off()


for (cancType in setdiff(unique(sample.ann[, 3]), c("LAML", "HNSC"))) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]
	
	vars = apply(complete.cna[, samples], 1, IQR, na.rm=TRUE)

	colv=as.dendrogram(hclust(as.dist(1-cor(complete.cna[names(sort(vars, decreasing=TRUE)[1:5000]), samples])), method="complete"))
	rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.cna[names(sort(vars, decreasing=TRUE)[1:5000]), samples]))), method="complete"))
	
	png(paste("tcga_ccla_cna_", cancType, ".png", sep=""), width=4000, height=3000, res=300)
	aheatmap(complete.cna[names(sort(vars, decreasing=TRUE)[1:5000]), samples],
         	 scale="none",
         	 color=colorpanel(49, low="blue", mid="white", high="red"),
         	 breaks=seq(min(complete.cna[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
         	 						max(complete.cna[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
         	 						length.out=50),
         	 Colv=colv,
         	 Rowv=rowv,
         	 distfun="pearson",
         	 annCol=sample.ann[samples, ],
         	 labRow=NA,
         	 labCol=NA)
  dev.off()


  vars = apply(complete.cna.combat[, samples], 1, IQR, na.rm=TRUE)

  colv=as.dendrogram(hclust(as.dist(1-cor(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples])), method="complete"))
  rowv=as.dendrogram(hclust(as.dist(1-cor(t(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples]))), method="complete"))
         
  png(paste("tcga_ccla_cna_", cancType, "_combat.png", sep=""), width=4000, height=3000, res=300)
  aheatmap(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples],
  	       scale="none",
   	       color=colorpanel(49, low="blue", mid="white", high="red"),
   	       breaks=seq(min(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						max(complete.cna.combat[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						length.out=50),
   	       Colv=colv,
   	       Rowv=rowv,
   	       distfun="pearson",
   	       annCol=sample.ann[samples, ],
   	       labRow=NA,
   	       labCol=NA)
  dev.off()
}

for (cancType in setdiff(unique(sample.ann[, 3]), c("LAML", "HNSC"))) {
	samples = rownames(sample.ann)[which(sample.ann[, 3] == cancType)]

	temp = complete.cna[, samples]
	mod = model.matrix(~as.factor(cancer), data=data.frame(sample.ann[samples, c("cancer", "tissue")]))								 
	# Normalize for TCGA vs CCLE
	temp = ComBat(temp, batch=sample.ann[, "batch"], mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

  vars = apply(temp, 1, IQR, na.rm=TRUE)

  colv=as.dendrogram(hclust(as.dist(1-cor(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples])), method="complete"))
  rowv=as.dendrogram(hclust(as.dist(1-cor(t(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples]))), method="complete"))
         
  png(paste("tcga_ccla_cna_", cancType, "_sepcombat.png", sep=""), width=4000, height=3000, res=300)
  aheatmap(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples],
  	       scale="none",
   	       color=colorpanel(49, low="blue", high="yellow"),
   	       breaks=seq(min(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						max(temp[names(sort(vars, decreasing=TRUE)[1:5000]), samples]), 
   	       						length.out=50),
   	       Colv=colv,
   	       Rowv=rowv,
   	       distfun="pearson",
   	       annCol=sample.ann[samples, ],
   	       labRow=NA,
   	       labCol=NA)
  dev.off()
}

