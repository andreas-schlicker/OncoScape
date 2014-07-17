# Data
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_exprs.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/DOWNLOAD_20120312/RDATA/cl_ccle_exprs.rdata")
# CCLE CL to TCGA mapping
cls.df = read.table("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/CCLE_TCGA_MAPPING/ccle_tcga_mapping.tsv", sep="\t", quote="", header=TRUE, stringsAsFactors=FALSE)
rownames(cls.df) = cls.df[, 2]
# Cell lines mapped to TCGA tissues that also have gene expression data
mappedCcle = intersect(as.character(cls.df[, "CL"]), colnames(cl.ccle.exprs))

# List of genes in all matrices
commonGenes = rownames(cl.ccle.exprs)
for (i in 1:length(tcga.exprs)) {
	commonGenes = intersect(commonGenes, rownames(tcga.exprs[[i]][[1]]))
}
for (i in 1:length(tcga.mn.exprs)) {
	commonGenes = intersect(commonGenes, rownames(tcga.mn.exprs[[i]][[1]]))
}

# Add a _mn to all matched normal expression samples
tcga.mn.exprs = lapply(tcga.mn.exprs, function(x) { colnames(x[[1]]) = paste(colnames(x[[1]]), "mn", sep="_"); x })


# Combine everything into one big expression matrix
complete.exprs = matrix(NA, nrow=length(commonGenes), ncol=length(mappedCcle)+sum(unlist(lapply(tcga.exprs, function(x) { ncol(x[[1]]) }))) + sum(unlist(lapply(tcga.mn.exprs, function(x) { ncol(x[[1]]) }))))
rownames(complete.exprs) = commonGenes
colnames(complete.exprs) = rep("", times=ncol(complete.exprs))

i = 1
last = i + length(mappedCcle) - 1
colnames(complete.exprs)[i:last] = mappedCcle
complete.exprs[, i:last] = cl.ccle.exprs[commonGenes, mappedCcle] - apply(cl.ccle.exprs[commonGenes, mappedCcle], 1, mean)

for (j in 1:length(tcga.exprs)) {
	i = last + 1
	last = i + ncol(tcga.exprs[[j]][[1]]) - 1
	colnames(complete.exprs)[i:last] = colnames(tcga.exprs[[j]][[1]])

	temp = tcga.exprs[[j]][[1]][commonGenes, ]
	
	if (min(temp, na.rm=TRUE) >= 0) {	
		temp = log2(temp + 1) - apply(log2(temp + 1), 1, mean)
	}

	complete.exprs[, i:last] = temp
}

for (j in 1:length(tcga.mn.exprs)) {
	i = last + 1
	last = i + ncol(tcga.mn.exprs[[j]][[1]]) - 1
	colnames(complete.exprs)[i:last] = colnames(tcga.mn.exprs[[j]][[1]])
	
	temp = tcga.mn.exprs[[j]][[1]][commonGenes, ]
	
	if (min(temp, na.rm=TRUE) >= 0) {	
		temp = log2(temp + 1) - apply(log2(temp + 1), 1, mean)
	}

	complete.exprs[, i:last] = temp	
}


# Sample annoation data.frame
sample.ann = rbind(data.frame(batch="ccle", cancer="cancer", tissue=cls.df[mappedCcle, 1]),
									 data.frame(batch="tcga", cancer="cancer", tissue=unlist(lapply(names(tcga.exprs), function(x) { rep(x, times=ncol(tcga.exprs[[x]][[1]])) }))),
									 data.frame(batch="tcga", cancer="normal", tissue=unlist(lapply(names(tcga.mn.exprs), function(x) { rep(x, times=ncol(tcga.mn.exprs[[x]][[1]])) }))))

# Replace the original matrices with the corrected ones
for (n in names(tcga.exprs)) {
	tcga.exprs[[n]][[1]] = complete.exprs[, colnames(tcga.exprs[[n]][[1]])]
}
for (n in names(tcga.mn.exprs)) {
	tcga.mn.exprs[[n]][[1]] = complete.exprs[, colnames(tcga.mn.exprs[[n]][[1]])]
	colnames(tcga.mn.exprs[[n]][[1]]) = gsub("_mn", "", colnames(tcga.mn.exprs[[n]][[1]])) 
}
# Create the same type of data structure for CCLE
ccle.exprs = list()
for (n in as.character(unique(cls.df[, "TCGA"]))) {
	ccle.exprs[[n]][[1]] = complete.exprs[, intersect(mappedCcle, as.character(cls.df[which(cls.df[, 1] == n), 2]))]
}

save(tcga.exprs, tcga.mn.exprs, ccle.exprs, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_exprs_ccle.rdata")

