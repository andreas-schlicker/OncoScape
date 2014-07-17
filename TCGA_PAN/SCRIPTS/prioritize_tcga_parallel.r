#!/usr/bin/Rscript --slave 

# Functions for loading TCGA data
source("/srv/nfs4/medoid-home/NKI/a.schlicker/R_SCRIPTS/tcga_synapse.r")

# Prioritization scripts
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/methylationmap.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/cnamap.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/utils.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/sommut.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/achilles.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/expression.r")
source("/srv/nfs4/medoid-home/NKI/a.schlicker/cagepR/scoring.r")

OPTIONS = getOptionList(c("config", "output", "genes", "genecol", "cancer", "samples", "threads", "methylation", "filtercol", "filterscore", "filtersnps"))

# Is there a configuration file?
args = commandArgs(TRUE)
if (length(args) == 0) {
	args = "--help"
}
config = parseOptions(OPTIONS, args)
if (!is.null(config$config)) {
	config = parseConfigFile(config$config)
} else {
	config = c()
}
# Parse options from configuration file and command line arguments
configOptions = parseOptions(OPTIONS, c(config, commandArgs(TRUE)))

# Read the gene list
genes = read.table(configOptions$genes, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)[, configOptions$genecol]

if (configOptions$samples != "") {
	SAMPLES = read.table(configOptions$samples, header=FALSE, sep="\t", quote="", as.is=TRUE)
} else {
	SAMPLES = NULL
}

if (length(configOptions$filtersnps) == 1 && configOptions$filtersnps == "none") {
	configOptions$filtersnps = NULL
}

if (!is.null(configOptions$filtercol) && configOptions$filtercol == "NULL") {
	configOptions$filtercol = NULL
}

if (length(configOptions$cancer) == 1 && configOptions$cancer == "all") {
	configOptions$cancer = setdiff(getCoreCancerTypes(), c("LAML"))
} else { 
	configOptions$cancer = intersect(configOptions$cancer, setdiff(getCoreCancerTypes(), c("LAML")))
}

configOptions$threads = min(configOptions$threads, length(configOptions$cancer))

library(foreach)
library(doMC)

registerDoMC()

# Load Achilles data
loadAchillesData("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/PROJECT_ACHILLES/20130620/RDATA")

# Load TCGA data
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_cna.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_exprs.rdata")
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_sommut.rdata")
# Load methylation data
if ("base" == configOptions$methylation) {
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_annotation.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_gene2probe.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_probe2gene.rdata")
} else if ("region" == configOptions$methylation) {
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth_region.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_annotation_region.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_gene2probe_region.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_probe2gene_region.rdata")
	configOptions$filtersnps = c()
} else {
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_meth_region_filtered.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_annotation_region_filtered.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_gene2probe_region_filtered.rdata")
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/METHYLATION/illumina_infinium450_probe2gene_region_filtered.rdata")
	configOptions$filtersnps = c()
}
filtered.probes450 = filterProbes(tcga.meth$LUSC, tcga.mn.meth$LUSC, genes, infinium450.probe.ann, gene2probe, snps=configOptions$filtersnps)
filtered.probes27 = filterProbes(tcga.meth$READ, tcga.mn.meth$READ, genes, infinium450.probe.ann, gene2probe, snps=configOptions$filtersnps)


type.mapping = list(LUSC=list(Type="Lung NSCLC", Subtype="Squamous"), 
				 					  READ=list(),
				 					  GBM=list(Type="GBM"),
				 					  BLCA=list(Type="Bladder"),
				 					  UCEC=list(Type="Endometrial"),
				 					  COAD=list(Type="Colon"),
				 					  OV=list(Type="Ovarian"),
				 					  LAML=list(Type="Leukemia"),
				 					  HNSC=list(),
				 					  LUAD=list(Type=c("Lung NSCLC", "Lung SCLC"), Subtype=c("Adenocarcinoma", "Adenocarcinoma, BAC")),
				 					  BRCA=list(Type="Breast"),
				 					  KIRC=list(Type="Renal Cell Carcinoma"))
				 					  
add.params = list(LUSC=list(platform="IlluminaHiSeq_RNASeqV2"),
									READ=list(platform="AgilentG4502A_07_3"),
									GBM=list(platform="AgilentG4502A_07_2"),
									BLCA=list(platform="IlluminaHiSeq_RNASeqV2"),
									UCEC=list(platform="IlluminaHiSeq_RNASeqV2"),
									COAD=list(platform="AgilentG4502A_07_3"),
									OV=list(platform="AgilentG4502A_07_3"),
									HNSC=list(platform="IlluminaHiSeq_RNASeqV2"),
									LUAD=list(platform="IlluminaHiSeq_RNASeqV2"),
									BRCA=list(platform="IlluminaHiSeq_RNASeqV2"),
									KIRC=list(platform="IlluminaHiSeq_RNASeqV2"))									

runCancer = function(cancType) {
	print(cancType)

	ct.exprs = tcga.exprs[[cancType]]
	ct.mn.exprs = tcga.mn.exprs[[cancType]]
	
	ct.cna = tcga.cna[[cancType]][[1]]
	ct.mn.cna = tcga.mn.cna[[cancType]][[1]]
	
	ct.meth = tcga.meth[[cancType]]
	ct.mn.meth = tcga.mn.meth[[cancType]]

	ct.sommut = tcga.sommut[[cancType]]

	#### Get the mapped cell lines
	CLS = unlist(lapply(type.mapping[[cancType]], function(x) { rownames(achilles.clann[which(achilles.clann[, "Type"] %in% x), ]) }))
	if (!is.null(type.mapping[[cancType]]$Subtype)) {
		CLS = intersect(CLS, rownames(achilles.clann[which(achilles.clann[, "Subtype"] %in% type.mapping[[cancType]]$Subtype), ]))
	}

	#### Check expression of the genes in normal tissue
	if (min(ct.exprs[[1]], na.rm=TRUE) >= 0) {	
		ct.exprs.logr = log2(ct.exprs[[1]] + 1) - apply(log2(ct.exprs[[1]] + 1), 1, mean)
		ct.mn.exprs.logr = log2(ct.mn.exprs[[1]] + 1) - apply(log2(ct.mn.exprs[[1]] + 1), 1, mean)
	} else {
		ct.exprs.logr = ct.exprs[[1]]
		ct.mn.exprs.logr = ct.mn.exprs[[1]]
	}

	exprs = doExprAnalysis(ct.exprs.logr, ct.mn.exprs.logr, genes, samples=SAMPLES, paired=TRUE)
	exprs.ts = summarizeExpr(ct.exprs.logr, ct.mn.exprs.logr, exprs, genes, wilcox.FDR=0.05, "down", paired=TRUE)
	exprs.og = summarizeExpr(ct.exprs.logr, ct.mn.exprs.logr, exprs, genes, wilcox.FDR=0.05, "up", paired=TRUE)
	
	#### Run Project Achilles analysis
	achilles.thresholds = quantile(achilles, probs=c(0.25, 0.75))
	
	achilles.og = doAchillesAnalysis(genes, score="phenotype", summarize=ltCutoffPercent(achilles.thresholds[1]), cls=CLS, relative=TRUE)
	achilles.og.sum = summarizeAchilles(achilles.og, 0.25, alternative="greater")
	
	achilles.ts = doAchillesAnalysis(genes, score="phenotype", summarize=gtCutoffPercent(achilles.thresholds[2]), cls=CLS, relative=TRUE)
	achilles.ts.sum = summarizeAchilles(achilles.ts, 0.25, alternative="greater")
	
	
	#### Run copy number analysis
	cgh = doCnaAnalysis(ct.cna, ct.mn.cna, ct.exprs.logr, genes, samples=SAMPLES)
	cgh.ts = summarizeCna(ct.cna, ct.mn.cna, cgh, genes, wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=-0.1,
	                      regulation="down", stddev=1)
	cgh.og = summarizeCna(ct.cna, ct.mn.cna, cgh, genes, wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=0.1,
	                      regulation="up", stddev=1)
	
  #### Run methylation analysis
	if (nrow(ct.meth) == 24980) {
		filtered.probes = filtered.probes27
	} else {
		filtered.probes = filtered.probes450
	}
	meth = doMethylationAnalysis(ct.meth, ct.mn.meth, ct.exprs.logr, probe2gene, filtered.probes, samples=SAMPLES)
	meth.ts = summarizeMethylation(ct.meth, ct.mn.meth, meth, gene2probe, probe2gene, filtered.probes, genes,
																 wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=0.1, regulation="down", gene.region=TRUE, stddev=1)
	meth.og = summarizeMethylation(ct.meth, ct.mn.meth, meth, gene2probe, probe2gene, filtered.probes, genes,
																 wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=-0.1, regulation="up", gene.region=TRUE, stddev=1)
	
	#### Perform mutation analysis
	ct.sommut[, "Tumor_Sample_Barcode"] = str_sub(ct.sommut[, "Tumor_Sample_Barcode"], start=1, end=12)
	sommut = doMutationAnalysis(ct.sommut, genes, ignore=c("Silent"), SAMPLES, filtercol=configOptions$filtercol, filter=configOptions$filterscore,
															genecol=1, typecol=9, samplecol=16, chromcol=5, 
															startcol=6, endcol=7, reference=11, tumor1=12, tumor2=13)
  sommut.ts = summarizeMutations(ct.sommut, sommut, genes, og.cutoff=0.2, ts.cutoff=0.2, ts.cutoff.low=0.05, 
							  								 score="ts", genecol=1, typecol=9, samplecol=16)
	sommut.og = summarizeMutations(ct.sommut, sommut, genes, og.cutoff=0.2, ts.cutoff=0.2, ts.cutoff.low=0.05, 
							  								 score="og", genecol=1, typecol=9, samplecol=16)
	
	#### Perform final prioritization
	## Complete list of tumor samples
	allTumors = sampleUnion(list(cna=colnames(ct.cna), 
															 exprs=colnames(ct.exprs),
															 meth=colnames(ct.meth),
															 sommut=ct.sommut[, "Tumor_Sample_Barcode"]))

	## Putative tumor suppressors
	prioritize.ts = cbind(ts.cna=cgh.ts$scores[genes],
												ts.methylation=meth.ts$scores[genes],
	                      ts.mutations=sommut.ts$scores[genes],
	                      ts.achilles=achilles.ts.sum[genes],
	                      ts.exprs=exprs.ts$scores[genes])
	
  ts.affected = sapply(genes, function(x) { length(sampleUnion(list(ts.cna=ifPresent(x, cgh.ts$samples),
																												 				    ts.methylation=ifPresent(x, meth.ts$samples),
																												 				    ts.mutations=ifPresent(x, sommut.ts$samples),
																												 				    ts.exprs=ifPresent(x, exprs.ts$samples)))) })

  prioritize.ts = cbind(ts.score=dataTypeScore(prioritize.ts)[rownames(prioritize.ts)], 
  											prioritize.ts,
  											ts.affected.abs=ts.affected[genes],
  											ts.affected.rel=ts.affected[genes] / length(allTumors))
	
	## Putative oncogenes
	prioritize.og = cbind(og.cna=cgh.og$scores[genes],
 	                      og.methylation=meth.og$scores[genes],
	                      og.mutations=sommut.og$scores[genes],
	                      og.achilles=achilles.og.sum[genes],
	                      og.exprs=exprs.og$scores[genes])
	
  og.affected = sapply(genes, function(x) { length(sampleUnion(list(og.cna=ifPresent(x, cgh.og$samples),
																												 				    og.methylation=ifPresent(x, meth.og$samples),
																												 				    og.mutations=ifPresent(x, sommut.og$samples),
																												 				    og.exprs=ifPresent(x, exprs.og$samples)))) })
	                           
  prioritize.og = cbind(og.score=dataTypeScore(prioritize.og)[rownames(prioritize.og)],
  											prioritize.og,
  											og.affected.abs=og.affected[genes],
  											og.affected.rel=og.affected[genes] / length(allTumors))
	
	# Calculate a combined score that penalizes genes with conflicting evidence
	combined.score = prioritize.og[rownames(prioritize.ts), "og.score"] - prioritize.ts[rownames(prioritize.ts), "ts.score"]
	
	# Combine all summaries and detailed scores
	prioritize.combined = cbind(combined.score=combined.score[rownames(prioritize.ts)],
	                            prioritize.og[rownames(prioritize.ts), ], prioritize.ts[rownames(prioritize.ts), ])
	prioritize.details = cbind(meth.smpl.abs.ts=sortAddMissing(genes, meth.ts$summary[, "absolute"]),
                             meth.smpl.rel.ts=sortAddMissing(genes, meth.ts$summary[, "relative"]),
                             meth.smpl.abs.og=sortAddMissing(genes, meth.og$summary[, "absolute"]),
                             meth.smpl.rel.og=sortAddMissing(genes, meth.og$summary[, "relative"]),
                             cgh.diff=sortAddMissing(genes, cgh.ts$cna.analysis$diffs),
                             cgh.diff.fdr=sortAddMissing(genes, cgh.ts$cna.analysis$wilcox[, "wilcox.FDR"]),
                             cgh.cor=sortAddMissing(genes, cgh.ts$cna.analysis$cors[, "cor"]),
                             cgh.cor.fdr=sortAddMissing(genes, cgh.ts$cna.analysis$cors[, "cor.FDR"]),
                             cgh.smpl.abs.ts=sortAddMissing(genes, cgh.ts$summary[, "absolute"]),
                             cgh.smpl.rel.ts=sortAddMissing(genes, cgh.ts$summary[, "relative"]),
                             cgh.smpl.abs.og=sortAddMissing(genes, cgh.og$summary[, "absolute"]),
                             cgh.smpl.rel.og=sortAddMissing(genes, cgh.og$summary[, "relative"]),
                             sommut.ts=sortAddMissing(genes, sommut[, "ts"]),
                             sommut.og=sortAddMissing(genes, sommut[, "og"]),
                             sommut.total.muts=sortAddMissing(genes, sommut[, "total.mutations"]),
                             sommut.unique.muts=sortAddMissing(genes, sommut[, "unique.mutations"]),
                             sommut.mutated.samples=sortAddMissing(genes, sommut[, "mutated.samples"]),
                             sommut.total.samples=sortAddMissing(genes, sommut[, "total.samples"]),
                             sommut.smpl.abs.ts=sortAddMissing(genes, sommut.ts$summary[, "absolute"]),
                             sommut.smpl.rel.ts=sortAddMissing(genes, sommut.ts$summary[, "relative"]),
                             sommut.smpl.abs.og=sortAddMissing(genes, sommut.og$summary[, "absolute"]),
                             sommut.smpl.rel.og=sortAddMissing(genes, sommut.og$summary[, "relative"]),
                             exprs.tumor=sortAddMissing(genes, exprs$exprs[, "tumor"]),
                             exprs.normal=sortAddMissing(genes, exprs$exprs[, "normal"]),
                             exprs.diff.fdr=sortAddMissing(genes, exprs.ts$expr.analysis$wilcox[, "wilcox.FDR"]),
                             exprs.smpl.abs.ts=sortAddMissing(genes, exprs.ts$summary[, "absolute"]),
                             exprs.smpl.rel.ts=sortAddMissing(genes, exprs.ts$summary[, "relative"]),
                             exprs.smpl.abs.og=sortAddMissing(genes, exprs.og$summary[, "absolute"]),
                             exprs.smpl.rel.og=sortAddMissing(genes, exprs.og$summary[, "relative"]))
	                                
  list(samples=SAMPLES, cls=CLS,
  		 exprs=exprs, exprs.ts=exprs.ts, exprs.og=exprs.og,
  		 achilles.og=achilles.og, achilles.og.sum=achilles.og.sum, achilles.ts=achilles.ts, achilles.ts.sum=achilles.ts.sum,
  		 cgh=cgh, cgh.ts=cgh.ts, cgh.og=cgh.og,
  		 meth=meth, meth.ts=meth.ts, meth.og=meth.og,
  		 sommut=sommut, sommut.ts=sommut.ts, sommut.og=sommut.og,
  		 allTumors=allTumors, ts.affected=ts.affected, og.affected=og.affected,
  		 prioritize.ts=prioritize.ts, prioritize.og=prioritize.og, combined.score=combined.score, prioritize.combined=prioritize.combined, prioritize.details=prioritize.details)
}

results = mclapply(configOptions$cancer, runCancer, mc.cores=configOptions$threads)
names(results) = configOptions$cancer

save(results, genes, type.mapping, file=configOptions$output)

