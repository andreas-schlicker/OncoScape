library(doMC)
library(foreach)

registerDoMC()

load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_annotation.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_gene2probe.rdata")

summarizeRegions = function(methMat, gene2probe) {
	lapply(gene2probe, function(x) { lapply(x, function(y) { y = intersect(y, rownames(methMat)) ; if (is.null(y) || is.na(y) || length(y) == 0) NA else apply(methMat[y, , drop=FALSE], 2, median, na.rm=TRUE) }) })
}

averageRegion = function(methData, gene2probe) {
	region.meth = summarizeRegions(methData, gene2probe)
	for (i in names(region.meth)) {
		for (j in names(region.meth[[i]])) {
			if (is.na(region.meth[[i]][[j]])) {
				region.meth[[i]][[j]] = NULL
			}
		}
	}
	
	samples = colnames(methData)
	methMat.regions = matrix(NA, nrow=sum(unlist(lapply(region.meth, length))), ncol=ncol(methData))
	colnames(methMat.regions) = colnames(methData)
	rn = character(nrow(methMat.regions))
	j = 1
	for (gene in names(region.meth)) {
		for (region in names(region.meth[[gene]])) {
			methMat.regions[j, samples] = region.meth[[gene]][[region]][samples]
			rn[j] = paste(gene, region, sep="_")
			j = j + 1
		}
	}
	rownames(methMat.regions) = rn
	
	methMat.regions
}

# Filter out all probes with SNPs
excl.snps = names(which(apply(infinium450.probe.ann[, c("SNP_target", "SNP_within_10", "SNP_outside_10", "SNP_probe"), drop=FALSE], 1, any)))
tcga.meth.filtered = mclapply(tcga.meth, function(x) { x[setdiff(rownames(x), excl.snps), , drop=FALSE] }, mc.cores=11)
tcga.mn.meth.filtered = mclapply(tcga.mn.meth, function(x) { x[setdiff(rownames(x), excl.snps), ,drop=FALSE] }, mc.cores=11)

# Summarize methylation for all samples excluding SNP probes
tcga.meth.region.filtered = mclapply(tcga.meth.filtered, averageRegion, gene2probe=gene2probe, mc.cores=11)
tcga.mn.meth.region.filtered = mclapply(tcga.mn.meth.filtered, averageRegion, gene2probe=gene2probe, mc.cores=11)

# Summarize methylation for all samples including SNP probes
tcga.meth.region = mclapply(tcga.meth, averageRegion, gene2probe=gene2probe, mc.cores=11)
tcga.mn.meth.region = mclapply(tcga.mn.meth, averageRegion, gene2probe=gene2probe, mc.cores=11)


tcga.meth = tcga.meth.region.filtered
tcga.mn.meth = tcga.mn.meth.region.filtered
save(tcga.meth, tcga.mn.meth, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth_region_filtered.rdata")

tcga.meth = tcga.meth.region
tcga.mn.meth = tcga.mn.meth.region
save(tcga.meth, tcga.mn.meth, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth_region.rdata")

