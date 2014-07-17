library(doMC)
library(foreach)

registerDoMC()

load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_annotation.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_gene2probe.rdata")

gene2probe.orig = gene2probe
gene2probe = lapply(gene2probe, function(x) { x$Body=NULL; unlist(x) })

summarizeRegions = function(methMat, gene2probe) {
	lapply(gene2probe, function(y) { y = intersect(y, rownames(methMat)) ; if (is.null(y) || is.na(y) || length(y) == 0) NA else apply(methMat[y, , drop=FALSE], 2, median, na.rm=TRUE) })
}

averageRegion = function(methData, gene2probe) {
	region.meth = summarizeRegions(methData, gene2probe)
	for (i in names(region.meth)) {
		if (is.na(region.meth[[i]])) {
				region.meth[[i]] = NULL
		}
	}
	
	samples = colnames(methData)
	methMat.regions = matrix(NA, nrow=length(region.meth), ncol=ncol(methData))
	rownames(methMat.regions) = names(region.meth)
	colnames(methMat.regions) = colnames(methData)
	for (gene in names(region.meth)) {
		methMat.regions[j, samples] = region.meth[[gene]][samples]
	}
	
	methMat.regions
}

# Summarize methylation for all samples including SNP probes
tcga.meth.gene = mclapply(tcga.meth, averageRegion, gene2probe=gene2probe, mc.cores=11)
tcga.mn.meth.gene = mclapply(tcga.mn.meth, averageRegion, gene2probe=gene2probe, mc.cores=11)


tcga.meth = tcga.meth.gene
tcga.mn.meth = tcga.mn.meth.gene
save(tcga.meth, tcga.mn.meth, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth_gene.rdata")

