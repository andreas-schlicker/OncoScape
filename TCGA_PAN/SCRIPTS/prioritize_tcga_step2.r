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

OPTIONS = getOptionList(c("config", "input", "inputstep2", "output", "genes", "genecol", "cancer", 
													"samples", "threads", "datatype", "methylation", "filtersnps", "cellline", "combat"))

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
if (configOptions$genes != "") {
	genes = read.table(configOptions$genes, header=TRUE, sep="\t", quote="", stringsAsFactors=FALSE)[, configOptions$genecol]
} else {
	genes = NULL
}

if (configOptions$samples != "") {
	SAMPLES = read.table(configOptions$samples, header=FALSE, sep="\t", quote="", as.is=TRUE)[, 1]
} else {
	SAMPLES = NULL
}

if (length(configOptions$cancer) == 1 && configOptions$cancer == "all") {
	configOptions$cancer = setdiff(getCoreCancerTypes(), c("LAML"))
} else { 
	configOptions$cancer = intersect(configOptions$cancer, setdiff(getCoreCancerTypes(), c("LAML")))
}

configOptions$threads = min(configOptions$threads, length(configOptions$cancer))

if (length(configOptions$datatype) == 1 && configOptions$datatype == "all") {
	configOptions$datatype = c("exprs", "achilles", "cgh", "meth", "sommut")
}

if (length(configOptions$filtersnps) == 1 && configOptions$filtersnps == "none") {
	configOptions$filtersnps = NULL
}

# Read the precomputed data
load(configOptions$input)
if (configOptions$inputstep2 != "") {
	load(configOptions$inputstep2)
}

combat = ifelse("TRUE" == toupper(configOptions$combat), TRUE, FALSE)
ccle = ifelse("ccle" == configOptions$cellline, TRUE, FALSE)

library(stringr)
library(foreach)
library(doMC)

registerDoMC()

# Load TCGA data
if ("cgh" %in% configOptions$datatype) {
	if (ccle) {
		# File with CCLE data and Combat corrected TCGA data
		load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_cna_ccle.rdata")
	} 
	if (!ccle || !combat) {
		load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_cna.rdata")
	}
	
	if (ccle) {
		if (combat) {
			# Set the Combat batch corrected objects
			tcga.cna = ccle.cna.combat
			tcga.mn.cna = tcga.mn.cna.combat
		} else {
			tcga.cna = ccle.cna
		}
	}
	
	# And delete all unnecessary objects 
	suppressWarnings(rm(tcga.cna.combat, tcga.mn.cna.combat, ccle.cna, ccle.cna.combat))
}
if ("meth" %in% configOptions$datatype) {
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
	
	filtered.probes27 = filterProbes(tcga.meth$READ, tcga.mn.meth$READ, genes=NULL, infinium450.probe.ann, gene2probe, snps=configOptions$filtersnps)
	filtered.probes450 = filterProbes(tcga.meth$LUSC, tcga.mn.meth$LUSC, genes=NULL, infinium450.probe.ann, gene2probe, snps=configOptions$filtersnps)
}
if ("sommut" %in% configOptions$datatype) {
	if ("ccle" == configOptions$cellline) {
		load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_sommut_ccle.rdata")
		tcga.sommut = ccle.sommut
	} else {
		load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_sommut.rdata")
	}
}
if ("achilles" %in% configOptions$datatype) {
	loadAchillesData("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/PROJECT_ACHILLES/20130620/RDATA")
}


# Load expression data
if (ccle) {
	# File with CCLE data and Combat corrected TCGA data
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_exprs_ccle.rdata")
} 
if (!ccle || !combat) {
	load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_exprs.rdata")
}
	
if (ccle) {
	if (combat) {
		# Set the Combat batch corrected objects
		tcga.exprs = ccle.exprs.combat
		tcga.mn.exprs = tcga.mn.exprs.combat
	} else {
		tcga.exprs = ccle.exprs
	}
}
	
# And delete all unnecessary objects 
suppressWarnings(rm(tcga.exprs.combat, tcga.mn.exprs.combat, ccle.exprs, ccle.exprs.combat))


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
	if (exists("tcga.exprs")) {
		ct.exprs = tcga.exprs[[cancType]]
		ct.mn.exprs = tcga.mn.exprs[[cancType]]
		
		if (min(ct.exprs[[1]], na.rm=TRUE) >= 0) {	
			ct.exprs.logr = log2(ct.exprs[[1]] + 1) - apply(log2(ct.exprs[[1]] + 1), 1, mean)
			ct.mn.exprs.logr = log2(ct.mn.exprs[[1]] + 1) - apply(log2(ct.mn.exprs[[1]] + 1), 1, mean)
		} else {
			ct.exprs.logr = ct.exprs[[1]]
			ct.mn.exprs.logr = ct.mn.exprs[[1]]
		}
	}
	
	if (exists("tcga.cna")) {
		ct.cna = tcga.cna[[cancType]][[1]]
		ct.mn.cna = tcga.mn.cna[[cancType]][[1]]
	}
	
	if (exists("tcga.meth")) {
		ct.meth = tcga.meth[[cancType]]
		ct.mn.meth = tcga.mn.meth[[cancType]]
	}

	if (exists("tcga.sommut")) {
		ct.sommut = tcga.sommut[[cancType]]
	}

	#### Check expression of the genes in normal tissue
	exprs = results.step1[[cancType]]$exprs
	if ("exprs" %in% configOptions$datatype) {
		exprs.ts = summarizeExpr(ct.exprs.logr, ct.mn.exprs.logr, exprs, genes, wilcox.FDR=0.05, "down", paired=TRUE)
		exprs.og = summarizeExpr(ct.exprs.logr, ct.mn.exprs.logr, exprs, genes, wilcox.FDR=0.05, "up", paired=TRUE)
	} else if (exists("results")) {
		exprs.ts = results[[cancType]]$exprs.ts
		exprs.og = results[[cancType]]$exprs.og
	} else {
		exprs.ts = NULL
		exprs.og = NULL
	}
	
	#### Run Project Achilles analysis
	achilles.thresholds = results.step1[[cancType]]$achilles.thresholds
	
	achilles.ts = results.step1[[cancType]]$achilles.ts
	achilles.og = results.step1[[cancType]]$achilles.og
	
	if ("achilles" %in% configOptions$datatype) {
		#### Get the mapped cell lines	
		CLS = unlist(lapply(type.mapping[[cancType]], function(x) { rownames(achilles.clann[which(achilles.clann[, "Type"] %in% x), ]) }))
		if (!is.null(type.mapping[[cancType]]$Subtype)) {
			CLS = intersect(CLS, rownames(achilles.clann[which(achilles.clann[, "Subtype"] %in% type.mapping[[cancType]]$Subtype), ]))
		}
		
		achilles.ts.sum = summarizeAchilles(achilles.ts, 0.25, alternative="greater")
		achilles.og.sum = summarizeAchilles(achilles.og, 0.25, alternative="greater")
	} else if (exists("results")) {
		CLS = results[[cancType]]$cls
		
		achilles.ts.sum = results[[cancType]]$achilles.ts.sum
		achilles.og.sum = results[[cancType]]$achilles.og.sum
	} else {
		CLS = NULL
		achilles.ts.sum = NULL
		achilles.og.sum = NULL
	}
	
	#### Run copy number analysis
	cgh = results.step1[[cancType]]$cgh
	if ("cgh" %in% configOptions$datatype) {
		cgh.ts = summarizeCna(ct.cna, ct.mn.cna, cgh, genes, wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=-0.1,
	  	                    regulation="down", stddev=1)
    cgh.og = summarizeCna(ct.cna, ct.mn.cna, cgh, genes, wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=0.1,
	  	                    regulation="up", stddev=1)
	} else if (exists("results")) {
		cgh.ts = results[[cancType]]$cgh.ts
		cgh.og = results[[cancType]]$cgh.og
	} else {
		cgh.ts = NULL
		cgh.og = NULL
	}
	
	#### Run methylation analysis
	meth = results.step1[[cancType]]$meth
	if ("meth" %in% configOptions$datatype) {
		if (nrow(ct.meth) < 30000) {
			filtered.probes = filtered.probes27
		} else {
			filtered.probes = filtered.probes450
		}

		meth.ts = summarizeMethylation(ct.meth, ct.mn.meth, meth, gene2probe.flat, probe2gene, filtered.probes, genes,
																	 wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=0.1, regulation="down", gene.region=TRUE, stddev=1)
	  meth.og = summarizeMethylation(ct.meth, ct.mn.meth, meth, gene2probe.flat, probe2gene, filtered.probes, genes,
																	 wilcox.FDR=0.05, cor.FDR=0.05, diff.cutoff=-0.1, regulation="up", gene.region=TRUE, stddev=1)
	} else if (exists("results")) {
		meth.ts = results[[cancType]]$meth.ts
		meth.og = results[[cancType]]$meth.og
	} else {
		meth.ts = NULL
		meth.og = NULL
	}
	
	#### Perform mutation analysis
	if ("sommut" %in% configOptions$datatype) {
		ct.sommut[, "Tumor_Sample_Barcode"] = str_sub(ct.sommut[, "Tumor_Sample_Barcode"], start=1, end=12)
		sommut = results.step1[[cancType]]$sommut
		sommut.ts = summarizeMutations(ct.sommut, sommut, genes, og.cutoff=0.2, ts.cutoff=0.2, ts.cutoff.low=0.05, 
								  								 score="ts", genecol=1, typecol=9, samplecol=16)
	  sommut.og = summarizeMutations(ct.sommut, sommut, genes, og.cutoff=0.2, ts.cutoff=0.2, ts.cutoff.low=0.05, 
								  								 score="og", genecol=1, typecol=9, samplecol=16)
	} else if (exists("results")) {
		sommut = results[[cancType]]$sommut
		sommut.ts = results[[cancType]]$sommut.ts
		sommut.og = results[[cancType]]$sommut.og
	} else {
		sommut = NULL
		sommut.ts = NULL
		sommut.og = NULL
	}
	
	#### Perform final prioritization
	## Complete list of tumor samples
	if (!is.null(SAMPLES)) {
		allTumors = SAMPLES
	} else {
		if (exists("results.step1") && !is.null(results.step1[[cancType]])) {
			allTumors = results.step1[[cancType]]$allTumors
		} else {
			allTumors = c()
		}
		allTumors = sampleUnion(list(allTumors=allTumors,
																 cna=if (exists("ct.cna")) colnames(ct.cna) else c(), 
							 			      		 	 exprs=if (exists("ct.exprs")) colnames(ct.exprs) else c(),
							 			      		 	 meth=if (exists("ct.meth")) colnames(ct.meth) else c(),
							 			      		 	 sommut=if (exists("ct.sommut")) ct.sommut[, "Tumor_Sample_Barcode"] else c()))
	}
	
  ## Putative tumor suppressors
	prioritize.ts = combineScores(list(ts.cna=cgh.ts$scores,
																		 ts.methylation=meth.ts$scores,
																		 ts.mutations=sommut.ts$scores,
																	   ts.achilles=achilles.ts.sum,
																		 ts.exprs=exprs.ts$scores),
																	genes)
		
  # Get all genes from the score matrix to make live easier later
  if (is.null(genes)) {
  	genes = rownames(prioritize.ts)
  }
  
  ts.affected = sapply(genes, function(x) { length(sampleUnion(list(ts.cna=ifPresent(x, cgh.ts$samples, allTumors),
																												 				    ts.methylation=ifPresent(x, meth.ts$samples, allTumors),
																												 				    ts.mutations=ifPresent(x, sommut.ts$samples, allTumors),
																												 				    ts.exprs=ifPresent(x, exprs.ts$samples, allTumors)))) })

  prioritize.ts = cbind(ts.score=dataTypeScore(prioritize.ts)[genes], 
  											prioritize.ts[genes, ],
  											ts.affected.abs=ts.affected[genes],
  											ts.affected.rel=ts.affected[genes] / length(allTumors))
	
	## Putative oncogenes
	prioritize.og = combineScores(list(og.cna=cgh.og$scores,
																		 og.methylation=meth.og$scores,
																		 og.mutations=sommut.og$scores,
																	   og.achilles=achilles.og.sum,
																		 og.exprs=exprs.og$scores),
																	genes)
																	
  og.affected = sapply(genes, function(x) { length(sampleUnion(list(og.cna=ifPresent(x, cgh.og$samples, allTumors),
																												 				    og.methylation=ifPresent(x, meth.og$samples, allTumors),
																												 				    og.mutations=ifPresent(x, sommut.og$samples, allTumors),
																												 				    og.exprs=ifPresent(x, exprs.og$samples, allTumors)))) })
	                           
  prioritize.og = cbind(og.score=dataTypeScore(prioritize.og)[genes],
  											prioritize.og[genes, ],
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
						 			      		 
	print(cancType)

	c(results.step1[[cancType]], list(samples=SAMPLES, cls=CLS, exprs.ts=exprs.ts, exprs.og=exprs.og, achilles.og.sum=achilles.og.sum, achilles.ts.sum=achilles.ts.sum, 
	  cgh.ts=cgh.ts, cgh.og=cgh.og, meth.ts=meth.ts, meth.og=meth.og, sommut.ts=sommut.ts, sommut.og=sommut.og, allTumors=allTumors, ts.affected=ts.affected, og.affected=og.affected,
	  prioritize.ts=prioritize.ts, prioritize.og=prioritize.og, combined.score=combined.score, prioritize.combined=prioritize.combined, prioritize.details=prioritize.details))
}

if (exists("results")) {
	results.loaded = results
}

results = mclapply(configOptions$cancer, runCancer, mc.cores=configOptions$threads)
names(results) = configOptions$cancer

if (exists("results.loaded")) {
	results = c(results, results.loaded[setdiff(names(results.loaded), names(results))])
}

save(results, file=configOptions$output)

