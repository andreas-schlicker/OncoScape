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

library(stringr)
library(foreach)
library(doMC)

registerDoMC()

OPTIONS = getOptionList(c("config", "input", "output", "cancer", "samples", "threads", "datatype", 
													"methylation", "filtercol", "filterscore", "filtersnps", "cellline", "combat"))

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
configOptions = parseOptions(OPTIONS, c(commandArgs(TRUE), config))

if (configOptions$samples != "") {
	SAMPLES = read.table(configOptions$samples, header=FALSE, sep="\t", quote="", as.is=TRUE)[, 1]
} else {
	SAMPLES = NULL
}

if (!is.null(configOptions$input) && configOptions$input != "") {
	load(configOptions$input)
} 

if (length(configOptions$cancer) == 1 && configOptions$cancer == "all") {
	configOptions$cancer = setdiff(getCoreCancerTypes(), c("LAML"))
} else { 
	configOptions$cancer = intersect(configOptions$cancer, setdiff(getCoreCancerTypes(), c("LAML")))
}

if (length(configOptions$datatype) == 1 && configOptions$datatype == "all") {
	configOptions$datatype = c("exprs", "achilles", "cgh", "meth", "sommut")
}

configOptions$threads = min(configOptions$threads, length(configOptions$cancer))

if (configOptions$filtercol == "NULL") {
	configOptions$filtercol = NULL
}

if (length(configOptions$filtersnps) == 1 && configOptions$filtersnps == "none") {
	configOptions$filtersnps = NULL
}

combat = ifelse("TRUE" == toupper(configOptions$combat), TRUE, FALSE)
ccle = ifelse("ccle" == configOptions$cellline, TRUE, FALSE)

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


genes = NULL

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
	if ("exprs" %in% configOptions$datatype) {
		exprs = doExprAnalysis(ct.exprs.logr, ct.mn.exprs.logr, genes, samples=SAMPLES, paired=TRUE)
	} else if (exists("results.step1")) {
		exprs = results.step1[[cancType]]$exprs
	} else {
		exprs = NULL
	}
	
	#### Run Project Achilles analysis
	if ("achilles" %in% configOptions$datatype) {
		#### Get the mapped cell lines
		CLS = unlist(lapply(type.mapping[[cancType]], function(x) { rownames(achilles.clann[which(achilles.clann[, "Type"] %in% x), ]) }))
		if (!is.null(type.mapping[[cancType]]$Subtype)) {
			CLS = intersect(CLS, rownames(achilles.clann[which(achilles.clann[, "Subtype"] %in% type.mapping[[cancType]]$Subtype), ]))
		}

		achilles.thresholds = quantile(achilles, probs=c(0.25, 0.75))
		if (!is.null(CLS)) {
			achilles.og = doAchillesAnalysis(genes, score="phenotype", summarize=ltCutoffPercent(achilles.thresholds[1]), cls=CLS)
			achilles.ts = doAchillesAnalysis(genes, score="phenotype", summarize=gtCutoffPercent(achilles.thresholds[2]), cls=CLS)
		} else {
			if (is.null(genes)) {
				g = rownames(achilles)
			}
			achilles.og = rep(0, times=length(g))
			names(achilles.og) = g
			achilles.ts = rep(0, times=length(g))
			names(achilles.ts) = g
		}
	} else if (exists("results.step1")) {
		CLS = results.step1[[cancType]]$cls
		achilles.thresholds = results.step1[[cancType]]$achilles.thresholds
		achilles.og = results.step1[[cancType]]$achilles.og
		achilles.ts = results.step1[[cancType]]$achilles.ts
	} else {
		CLS = NULL
		achilles.thresholds = NULL
		achilles.og = NULL
		achilles.ts = NULL
	}
	
	#### Run copy number analysis
	if ("cgh" %in% configOptions$datatype) {
		cgh = doCnaAnalysis(ct.cna, ct.mn.cna, ct.exprs.logr, genes, samples=SAMPLES)
	} else if (exists("results.step1")) {
		cgh = results.step1[[cancType]]$cgh
	} else {
		cgh = NULL
	}
	
	#### Run methylation analysis
	if ("meth" %in% configOptions$datatype) {
		if (nrow(ct.meth) == 24980) {
			filtered.probes = filtered.probes27
		} else {
			filtered.probes = filtered.probes450
		}
	
		meth = doMethylationAnalysis(ct.meth, ct.mn.meth, ct.exprs.logr, probe2gene.flat, filtered.probes, samples=SAMPLES)
	} else if (exists("results.step1")) {
		meth = results.step1[[cancType]]$meth
	} else {
		meth = NULL
	}
	
	#### Perform mutation analysis
	if ("sommut" %in% configOptions$datatype) {
		ct.sommut[, "Tumor_Sample_Barcode"] = str_sub(ct.sommut[, "Tumor_Sample_Barcode"], start=1, end=12)
		sommut = doMutationAnalysis(ct.sommut, genes, ignore=c("Silent"), SAMPLES, filtercol=configOptions$filtercol, filter=configOptions$filterscore,
																genecol=1, typecol=9, samplecol=16, chromcol=5, 
																startcol=6, endcol=7, reference=11, tumor1=12, tumor2=13)
	} else if (exists("results.step1")) {
		sommut = results.step1[[cancType]]$sommut
	} else {
		sommut = NULL
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
							 			      		 
	print(cancType)

	list(samples=SAMPLES, cls=CLS, exprs=exprs, achilles.thresholds=achilles.thresholds, achilles.og=achilles.og, achilles.ts= achilles.ts, 
	     cgh=cgh, meth=meth, sommut=sommut, allTumors=allTumors)
}

if (exists("results.step1")) {
	results.step1.loaded = results.step1
}

results.step1 = mclapply(configOptions$cancer, runCancer, mc.cores=configOptions$threads)
names(results.step1) = configOptions$cancer

if (exists("results.step1.loaded")) {
	results.step1 = c(results.step1, results.step1.loaded[setdiff(names(results.step1.loaded), names(results.step1))])
}

save(results.step1, file=configOptions$output)

