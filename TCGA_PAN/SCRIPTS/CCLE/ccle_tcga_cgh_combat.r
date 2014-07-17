# ComBat library
library(sva)
library(impute)

# Data
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_cna.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/DOWNLOAD_20120312/RDATA/cl_ccle_cna.rdata")
# CCLE CL to TCGA mapping
cls.df = read.table("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/CCLE_TCGA_MAPPING/ccle_tcga_mapping.tsv", sep="\t", quote="", header=TRUE, stringsAsFactors=FALSE)
rownames(cls.df) = cls.df[, 2]

# Cell lines mapped to TCGA tissues that also have gene expression data
mappedCcle = intersect(as.character(cls.df[, "CL"]), colnames(cl.ccle.cna))
# List of genes in all matrices
commonGenes = rownames(cl.ccle.cna)
for (i in 1:length(tcga.cna)) {
	commonGenes = intersect(commonGenes, rownames(tcga.cna[[i]][[1]]))
}
for (i in 1:length(tcga.mn.cna)) {
	commonGenes = intersect(commonGenes, rownames(tcga.mn.cna[[i]][[1]]))
}

# Add a _mn to all matched normal expression samples
tcga.mn.cna = lapply(tcga.mn.cna, function(x) { colnames(x[[1]]) = paste(colnames(x[[1]]), "mn", sep="_"); x })

ccle.cna = list()
ccle.cna.combat = list()
tcga.cna.combat = list()
tcga.mn.cna.combat = list()

# Combine everything into one big expression matrix
complete.cna = matrix(NA, nrow=length(commonGenes), ncol=length(mappedCcle)+sum(unlist(lapply(tcga.cna, function(x) { ncol(x[[1]]) }))) + sum(unlist(lapply(tcga.mn.cna, function(x) { ncol(x[[1]]) }))))
rownames(complete.cna) = commonGenes
colnames(complete.cna) = rep("", times=ncol(complete.cna))

i = 1
last = i + length(mappedCcle) - 1
colnames(complete.cna)[i:last] = mappedCcle
complete.cna[, i:last] = cl.ccle.cna[commonGenes, mappedCcle]

for (j in 1:length(tcga.cna)) {
	i = last + 1
	last = i + ncol(tcga.cna[[j]][[1]]) - 1
	colnames(complete.cna)[i:last] = colnames(tcga.cna[[j]][[1]])

	temp = tcga.cna[[j]][[1]][commonGenes, ]
	if (length(which(is.na(temp))) > 0) {
		temp = impute.knn(temp)$data
	}
	
	if (min(temp, na.rm=TRUE) >= 0) {	
		temp = log2(temp + 1) - apply(log2(temp + 1), 1, mean)
	}

	complete.cna[, i:last] = temp
}

for (j in 1:length(tcga.mn.cna)) {
	i = last + 1
	last = i + ncol(tcga.mn.cna[[j]][[1]]) - 1
	colnames(complete.cna)[i:last] = colnames(tcga.mn.cna[[j]][[1]])
	
	temp = tcga.mn.cna[[j]][[1]][commonGenes, ]
	if (length(which(is.na(temp))) > 0) {
		temp = impute.knn(temp)$data
	}
	
	if (min(temp, na.rm=TRUE) >= 0) {	
		temp = log2(temp + 1) - apply(log2(temp + 1), 1, mean)
	}

	complete.cna[, i:last] = temp
}


# Sample annoation data.frame
sample.ann = rbind(data.frame(batch="ccle", cancer="cancer", tissue=cls.df[mappedCcle, 1]),
									 data.frame(batch="tcga", cancer="cancer", tissue=unlist(lapply(names(tcga.cna), function(x) { rep(x, times=ncol(tcga.cna[[x]][[1]])) }))),
									 data.frame(batch="tcga", cancer="normal", tissue=unlist(lapply(names(tcga.mn.cna), function(x) { rep(x, times=ncol(tcga.mn.cna[[x]][[1]])) }))))
# Design matrix keeping cancer status and tissue type
mod = model.matrix(~as.factor(cancer)+as.factor(tissue), data=data.frame(sample.ann[, c("cancer", "tissue")]))								 
# Normalize for TCGA vs CCLE
complete.cna.combat = ComBat(complete.cna, batch=sample.ann[, "batch"], mod=mod, numCovs=NULL, par.prior=TRUE, prior.plots=FALSE)

# Replace the original matrices with the corrected ones
for (n in names(tcga.cna)) {
	tcga.cna.combat[[n]][[1]] = complete.cna.combat[, colnames(tcga.cna[[n]][[1]])]
	
	tcga.mn.cna.combat[[n]][[1]] = complete.cna.combat[, colnames(tcga.mn.cna[[n]][[1]])]
	colnames(tcga.mn.cna.combat[[n]][[1]]) = gsub("_mn", "", colnames(tcga.mn.cna.combat[[n]][[1]])) 

	ccle.cna[[n]][[1]] = complete.cna[, intersect(mappedCcle, as.character(cls.df[which(cls.df[, 1] == n), 2]))]
	ccle.cna.combat[[n]][[1]] = complete.cna.combat[, intersect(mappedCcle, as.character(cls.df[which(cls.df[, 1] == n), 2]))]
}

save(tcga.cna.combat, tcga.mn.cna.combat, ccle.cna, ccle.cna.combat, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_cna_ccle.rdata")

