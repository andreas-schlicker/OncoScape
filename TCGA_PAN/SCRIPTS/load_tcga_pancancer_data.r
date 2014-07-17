source("/srv/nfs4/medoid-home/NKI/a.schlicker/R_SCRIPTS/tcga_synapse.r")

# Load library
library(synapseClient)
library(stringr)
library(ggplot2)

# and login
synapseLogin()

# Set the cache directory
synapseCacheDir("/srv/nfs4/void-bulk/NKI/a.schlicker/TCGA/SYNAPSE")
				 					  
add.params = list(LUSC=list(exprs=list(platform="IlluminaHiSeq_RNASeqV2"), meth=list(platform="HumanMethylation450")),
									READ=list(exprs=list(platform="AgilentG4502A_07_3"), meth=list(platform="HumanMethylation450")),
									GBM=list(exprs=list(platform="AgilentG4502A_07_2"), meth=list(platform="HumanMethylation450")),
									BLCA=list(exprs=list(platform="IlluminaHiSeq_RNASeqV2"), meth=list(platform="HumanMethylation450")),
									UCEC=list(exprs=list(platform="IlluminaHiSeq_RNASeqV2"), meth=list(platform="HumanMethylation450")),
									COAD=list(exprs=list(platform="AgilentG4502A_07_3"), meth=list(platform="HumanMethylation27")),
									OV=list(exprs=list(platform="AgilentG4502A_07_3"), meth=list(platform="HumanMethylation27")),
									HNSC=list(exprs=list(platform="IlluminaHiSeq_RNASeqV2"), meth=list(platform="HumanMethylation450")),
									LUAD=list(exprs=list(platform="IlluminaHiSeq_RNASeqV2"), meth=list(platform="HumanMethylation450")),
									BRCA=list(exprs=list(platform="IlluminaHiSeq_RNASeqV2"), meth=list(platform="HumanMethylation450")),
									KIRC=list(exprs=list(platform="IlluminaHiSeq_RNASeqV2"), meth=list(platform="HumanMethylation450")))

tcga.exprs = list()
tcga.mn.exprs = list()
tcga.cna = list()
tcga.mn.cna = list()
tcga.meth = list()
tcga.mn.meth = list()
tcga.sommut = list()

for (cancType in setdiff(getCoreCancerTypes(), c("LAML"))) {
	print(cancType)
	ct.synid = getSynId(getCancerParams(cancType))
	
	rowN = ifelse(add.params[[cancType]]$exprs == "IlluminaHiSeq_RNASeqV2", "split", "keep")
	tcga.exprs[[cancType]] = getMyData(ct.synid, "geneexp", "tumor", rowN, addParams=add.params[[cancType]]$exprs)
	tcga.mn.exprs[[cancType]] = getMyData(ct.synid, "geneexp", "normal", rowN, addParams=add.params[[cancType]]$exprs)
}

for (cancType in setdiff(getCoreCancerTypes(), c("LAML"))) {
	print(cancType)
	ct.synid = getSynId(getCancerParams(cancType))
	
	tcga.cna[[cancType]] = getMyData(ct.synid, "cna", "tumor")
	tcga.mn.cna[[cancType]] = getMyData(ct.synid, "cna", "normal")
}

for (cancType in setdiff(getCoreCancerTypes(), c("LAML"))) {
for (cancType in c("GBM", "OV", "LUAD", "BRCA", "READ")) {
	print(cancType)
	ct.synid = getSynId(getCancerParams(cancType))

	tcga.meth[[cancType]] = getMyData(ct.synid, "methylation", "tumor", rowN="keep", addParams=add.params[[cancType]]$meth)
	tcga.mn.meth[[cancType]] = getMyData(ct.synid, "methylation", "normal", rowN="keep", addParams=add.params[[cancType]]$meth)
	#tcga.meth[[cancType]] = getMyData(ct.synid, "methylation", "tumor", rowN="keep", addParams=add.params[[cancType]])
	#tcga.mn.meth[[cancType]] = getMyData(ct.synid, "methylation", "normal", rowN="keep", addParams=add.params[[cancType]])
}

for (cancType in setdiff(getCoreCancerTypes(), c("LAML"))) {
	print(cancType)
	ct.synid = getSynId(getCancerParams(cancType))
	
	tcga.sommut[[cancType]] = getMyMutationData(cancType)
}

synapseLogout()

# Separate COAD and READ mutations
all.coad = unique(c(str_sub(colnames(tcga.cna$COAD[[1]]), start=1, end=12), 
										str_sub(colnames(tcga.exprs$COAD[[1]]), start=1, end=12), 
										str_sub(colnames(tcga.meth$COAD), start=1, end=12)))
all.read = unique(c(str_sub(colnames(tcga.cna$READ[[1]]), start=1, end=12), 
									  str_sub(colnames(tcga.exprs$READ[[1]]), start=1, end=12), 
									  str_sub(colnames(tcga.meth$READ), start=1, end=12)))

tcga.sommut$COAD = subset(tcga.sommut$COAD, str_sub(Tumor_Sample_Barcode, start=1, end=12) %in% all.coad)
tcga.sommut$READ = subset(tcga.sommut$READ, str_sub(Tumor_Sample_Barcode, start=1, end=12) %in% all.read)
										

ds.size = matrix(0, nrow=9, ncol=length(tcga.exprs))
colnames(ds.size) = names(tcga.exprs)
rownames(ds.size) = c("Exprs_t", "Exprs_n", "CNV_t", "CNV_n", "Meth_t", "Meth_n", "Mut", "Exprs_CNV_t", "Exprs_Meth_t")
ds.size.df = data.frame()
for (cancType in names(tcga.exprs)) {
	ds.size["Exprs_t", cancType] = ncol(tcga.exprs[[cancType]][[1]])
	ds.size["Exprs_n", cancType] = ncol(tcga.mn.exprs[[cancType]][[1]])
	
	ds.size["CNV_t", cancType] = ncol(tcga.cna[[cancType]][[1]])
	ds.size["CNV_n", cancType] = ncol(tcga.mn.cna[[cancType]][[1]])
	
	ds.size["Meth_t", cancType] = ncol(tcga.meth[[cancType]])
	ds.size["Meth_n", cancType] = ncol(tcga.mn.meth[[cancType]])
	
	ds.size["Mut", cancType] = length(unique(str_sub(tcga.sommut[[cancType]][, "Tumor_Sample_Barcode"], start=1, end=12)))
	
	ds.size["Exprs_CNV_t", cancType] = length(intersect(colnames(tcga.exprs[[cancType]][[1]]), colnames(tcga.cna[[cancType]][[1]])))
	ds.size["Exprs_Meth_t", cancType] = length(intersect(colnames(tcga.exprs[[cancType]][[1]]), colnames(tcga.meth[[cancType]])))
	
	ds.size.df = rbind(ds.size.df,
										 data.frame(cancer=cancType,
										 						dataset="Gene expression",
										 						samples=c(ds.size["Exprs_t", cancType], ds.size["Exprs_n", cancType]),
										 						sampletype=c("Tumor", "Normal")),
										 data.frame(cancer=cancType,
										 						dataset="Copy number",
										 						samples=c(ds.size["CNV_t", cancType], ds.size["CNV_n", cancType]),
										 						sampletype=c("Tumor", "Normal")),
										 data.frame(cancer=cancType,
										 						dataset="DNA methylation",
										 						samples=c(ds.size["Meth_t", cancType], ds.size["Meth_n", cancType]),
										 						sampletype=c("Tumor", "Normal")),
										 data.frame(cancer=cancType,
										 						dataset="Sommatic mutations",
										 						samples=c(ds.size["Mut", cancType], 0),
										 						sampletype=c("Tumor", "Normal")),
										 data.frame(cancer=cancType,
										 						dataset="Gene expression/Copy number",
										 						samples=ds.size["Exprs_CNV_t", cancType],
										 						sampletype="Tumor"),
										 data.frame(cancer=cancType,
										 						dataset="Gene expression/DNA methylation",
										 						samples=ds.size["Exprs_Meth_t", cancType],
										 						sampletype="Tumor"))
}

# Select methylation data sets with the largest overlap 
tcga.meth.all = tcga.meth
tcga.mn.meth.all = tcga.mn.meth

tcga.meth$LUSC = tcga.meth$LUSC[[2]]
tcga.mn.meth$LUSC = tcga.mn.meth$LUSC[[2]]
tcga.meth$READ = tcga.meth$READ[[1]]
tcga.mn.meth$READ = tcga.mn.meth$READ[[1]]
tcga.meth$GBM = tcga.meth$GBM[[2]]
tcga.mn.meth$GBM = tcga.mn.meth$GBM[[1]]
tcga.meth$BLCA = tcga.meth$BLCA[[1]]
tcga.mn.meth$BLCA = tcga.mn.meth$BLCA[[1]]
tcga.meth$UCEC = tcga.meth$UCEC[[2]]
tcga.mn.meth$UCEC = tcga.mn.meth$UCEC[[2]]
tcga.meth$COAD = tcga.meth$COAD[[1]]
tcga.mn.meth$COAD = tcga.mn.meth$COAD[[1]]
tcga.meth$OV = tcga.meth$OV[[1]]
tcga.mn.meth$OV = tcga.mn.meth$OV[[1]]
tcga.meth$HNSC = tcga.meth$HNSC[[1]]
tcga.mn.meth$HNSC = tcga.mn.meth$HNSC[[1]]
tcga.meth$LUAD = tcga.meth$LUAD[[2]]
tcga.mn.meth$LUAD = tcga.mn.meth$LUAD[[2]]
tcga.meth$BRCA = tcga.meth$BRCA[[2]]
tcga.mn.meth$BRCA = tcga.mn.meth$BRCA[[2]]
tcga.meth$KIRC = tcga.meth$KIRC[[2]]
tcga.mn.meth$KIRC = tcga.mn.meth$KIRC[[2]]


save(tcga.exprs, tcga.mn.exprs, file="tcga_pancancer4_exprs.rdata")
save(tcga.cna, tcga.mn.cna, file="tcga_pancancer4_cna.rdata")
save(tcga.meth, tcga.mn.meth, file="tcga_pancancer4_meth.rdata")
save(tcga.meth.all, tcga.mn.meth.all, file="tcga_pancancer4_meth_all.rdata")
save(tcga.sommut, file="tcga_pancancer4_sommut.rdata")
save(ds.size, ds.size.df, file="tcga_pancancer4_dssize.rdata")


png("dataset_sizes_pancancer.png", width=3000, height=3000, res=300)
ggplot(ds.size.df, aes(x=dataset, y=samples, fill=sampletype)) +
geom_bar(position=position_dodge(), stat="identity") +
coord_flip() +
facet_grid(cancer ~ sampletype) + 
ylab("Number of samples") + 
xlab("")
dev.off()


cols = c("#E69F00", "#888888", "#56B4E9", "#009E73")
dst = c("Gene expression", "Copy number", "DNA methylation", "Sommatic mutations")
names(dst) = c("Expr", "CNA", "Meth", "Mut")

# Sort the cancer types
ds.size.df$cancer = factor(ds.size.df$cancer, levels=sort(levels(ds.size.df$cancer)))

for (n in 1:length(dst)) {
	png(paste("dataset_sizes_pancancer_", tolower(names(dst)[n]), ".png", sep=""), width=4000, height=2000, res=300)
	
	p = ggplot(subset(ds.size.df, dataset==dst[n]), aes(x=cancer, y=samples, fill=dataset)) +
	geom_bar(position=position_dodge(), stat="identity") +
	scale_fill_manual(values=c(cols[n], "white")) +
	guides(fill=FALSE) +
	facet_grid(. ~ sampletype) +  
	ylab("Number of samples") + 
	xlab("Cancer Type") + 
	theme(axis.text.x=element_text(color="grey20", size=24, face="bold", angle=45, vjust=0.5),
				axis.text.y=element_text(color="grey20", size=24, face="bold"),
				axis.title.x=element_text(color="grey20", size=24, face="bold", vjust=0),
				axis.title.y=element_text(color="grey20", size=24, face="bold", vjust=0.3),
				strip.text=element_text(color="grey20", size=24, face="bold"))
	
	if (n %in% 2:3) {
		p = p + geom_point(data=subset(ds.size.df, dataset==paste("Gene expression/", dst[n], sep="")), 
									 		 aes(x=cancer, y=samples), color="white", size=6)
	}
	
	print(p)
	dev.off()
}

