# Load libraries 
library(biomaRt)
library(gplots)
library(NMF)

# Prioritization results
load("~/bulk/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140302/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata")

# Genes in all cancer types
genes = rownames(results$LUSC$prioritize.combined)
for (n in setdiff(names(results), "LUSC")) {
  genes = intersect(genes, rownames(results[[n]]$prioritize.combined))
}

# Combined the prioritization results across cancers
results.combined = results$LUSC$prioritize.combined[genes, ]
results.combined[which(is.na(results.combined))] = 0
for (n in setdiff(names(results), "LUSC")) {
	temp = results[[n]]$prioritize.combined[genes, ]
	temp[which(is.na(temp))] = 0
	results.combined = results.combined + temp
}
results.combined = results.combined / length(results)

# Find all gene locations
mart = useMart(host="ensembl.org", path="/biomart/martservice", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
geneLoc = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position", "strand", "band"), filters=c("hgnc_symbol"), values=genes, mart=mart)
# Remove all non-standard chromosome names and make everything numeric
geneLoc = geneLoc[which(geneLoc[, "chromosome_name"] %in% c("1", "2", "3", "4", "5", "6", "7", "8", 
																														"9", "10", "11", "12", "13", "14", "15", 
																														"16", "17", "18", "19", "20", "21", "22", 
																														"X", "Y")), ]
geneLoc[which(geneLoc[, "chromosome_name"] == "X"), "chromosome_name"] = "23"
geneLoc[which(geneLoc[, "chromosome_name"] == "Y"), "chromosome_name"] = "24"
geneLoc[, "chromosome_name"] = as.integer(geneLoc[, "chromosome_name"])

geneLoc[which(geneLoc[, "strand"] == 1), "strand"] = "+"
geneLoc[which(geneLoc[, "strand"] == -1), "strand"] = "-"

# Sort the gene information
geneLoc = geneLoc[order(geneLoc$chromosome_name, geneLoc$start_position), ]

chrom2genes = split(geneLoc$hgnc_symbol, geneLoc$chromosome_name)
names(chrom2genes) = paste("chr", names(chrom2genes), sep="")

medianScores = matrix(NA, ncol=3, nrow=24)
colnames(medianScores) = c("TS", "Combined", "OG")
rownames(medianScores) = paste("chr", 1:24, sep="")
medianScores[, "Combined"] = unlist(lapply(chrom2genes, function(x) { median(results.combined[x, "combined.score"], na.rm=TRUE) }))
medianScores[, "OG"] = unlist(lapply(chrom2genes, function(x) { median(results.combined[x, "og.score"], na.rm=TRUE) }))
medianScores[, "TS"] = -1*unlist(lapply(chrom2genes, function(x) { median(results.combined[x, "ts.score"], na.rm=TRUE) }))

png("median_scores_per_chromosome.png", width=3000, height=3000, res=300)
aheatmap(medianScores, scale="none", Colv=NA, Rowv=NA,
				 col=colorpanel(20, low="#034b87", mid="white", high="#880000"), breaks=0,
				 fontsize=20, fontface="bold")
dev.off()


meanScores = matrix(NA, ncol=3, nrow=24)
colnames(meanScores) = c("TS", "Combined", "OG")
rownames(meanScores) = paste("chr", 1:24, sep="")
meanScores[, "Combined"] = unlist(lapply(chrom2genes, function(x) { mean(results.combined[x, "combined.score"], na.rm=TRUE) }))
meanScores[, "OG"] = unlist(lapply(chrom2genes, function(x) { mean(results.combined[x, "og.score"], na.rm=TRUE) }))
meanScores[, "TS"] = -1*unlist(lapply(chrom2genes, function(x) { mean(results.combined[x, "ts.score"], na.rm=TRUE) }))

png("mean_scores_per_chromosome.png", width=3000, height=3000, res=300)
aheatmap(meanScores, scale="none", Colv=NA, Rowv=NA,
				 col=colorpanel(20, low="#034b87", mid="white", high="#880000"), breaks=0,
				 fontsize=20, fontface="bold")
dev.off()

