load("cl_ccle_muts.rdata")

cl.ccle.muts2 = cbind(cl.ccle.muts, cl.ccle.muts[, 9])
colnames(cl.ccle.muts2)[ncol(cl.ccle.muts2)] = "Variant_Classification_orig"
cl.ccle.muts2[which(cl.ccle.muts[, 9] == "Stop_Codon_Ins"), 9] = "Nonsense_Mutation"
cl.ccle.muts2[which(cl.ccle.muts[, 9] == "Splice_Site_Del"), 9] = "Splice_Site"
cl.ccle.muts2[which(cl.ccle.muts[, 9] == "Splice_Site_DNP"), 9] = "Splice_Site"
cl.ccle.muts2[which(cl.ccle.muts[, 9] == "Splice_Site_Ins"), 9] = "Splice_Site"
cl.ccle.muts2[which(cl.ccle.muts[, 9] == "Splice_Site_SNP"), 9] = "Splice_Site"

cls.df = read.table("/srv/nfs4/medoid-bulk/NKI/a.schlicker/CL/CCLE/CCLE_TCGA_MAPPING/ccle_tcga_mapping.tsv", sep="\t", quote="", header=TRUE, stringsAsFactors=FALSE)
rownames(cls.df) = cls.df[, 2]
mappedCcle = intersect(as.character(cls.df[, "CL"]), cl.ccle.muts2[, "Tumor_Sample_Barcode"])

ccle.sommut = list()
for (n in as.character(unique(cls.df[, "TCGA"]))) {
	ccle.sommut[[n]] = cl.ccle.muts2[which(cl.ccle.muts2[, "Tumor_Sample_Barcode"] %in% rownames(cls.df[which(cls.df[, "TCGA"] == n), ])), ]
}

save(ccle.sommut, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_sommut_ccle.rdata")

