setwd("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/SCORE_HEATMAPS")

# Load libraries 
library(biomaRt)
library(gplots)
library(NMF)

# Prioritization results
load("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata")

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

# Sort the genes according to location
sortedGenes = geneLoc[, "hgnc_symbol"]

## NOT RUN
## Manual processing to transform the heatmap into contour lines
## GIMP: Edge detect
# Plot the contour lines for the chromosomes
#sortedChrom = geneLoc[order(geneLoc$chromosome_name, geneLoc$start_position), ][, "chromosome_name"]
#sortedChromMat = matrix(c(sortedChrom, rep(0, times=(139^2)-length(sortedChrom))), nrow=139, byrow=FALSE)
#png("chrom_contours.png", width=4000, height=3000, res=300)
#aheatmap(sortedChromMat[, 1:138],
#				 color=colorpanel(24, low="grey80", high="blue", mid="darkgreen"),
#				 breaks=seq(1, 24, length.out=25),
#				 scale="none",
#				 Rowv=NA, Colv=NA, 
#				 cexRow=0, cexCol=0)
#dev.off()


## NOT RUN
## Manual processing to transform the heatmap into contour lines
## GIMP: Edge detect
# Plot the contour lines for the chromosome arms
#sortedBand = geneLoc[order(geneLoc$chromosome_name, geneLoc$start_position), ][, "band"]
#sortedBand[sortedBand == "p"] = 1
#sortedBand[sortedBand == "q"] = 2
#sortedBand = as.integer(sortedBand)
#sortedBandMat = matrix(c(sortedBand, rep(0, times=(139^2)-length(sortedBand))), nrow=139, byrow=FALSE)
#png("band_contours.png", width=4000, height=3000, res=300)
#aheatmap(sortedBandMat[, 1:138],
#				 color="black",
#				 scale="none",
#				 Rowv=NA, Colv=NA, 
#				 cexRow=0, cexCol=0)
#dev.off()


## Plot the different score heatmaps and add the chromosome contour lines
# For the averaged scores
for (n in setdiff(colnames(results.combined), c("og.affected.abs", "og.affected.rel", "ts.affected.abs", "ts.affected.rel"))) {
	sortedScores = matrix(c(results.combined[sortedGenes, n], rep(0, times=(139^2)-length(sortedGenes))), nrow=139, byrow=FALSE)[, 1:138]
	if (n == "combined.score") {
		cols = colorpanel(29, low="#034b87", mid="white", high="#880000")
		ext = max(abs(min(sortedScores, na.rm=TRUE)), max(sortedScores, na.rm=TRUE))
		#brks=seq(-1*ext, ext, length.out=30)
		brks = 0
	} else if (length(grep("ts.", n)) > 0) {
		cols = colorpanel(29, low="white", high="#034b87")
		brks=seq(0, max(sortedScores, na.rm=TRUE), length.out=30)
	} else {
		cols = colorpanel(29, low="white", high="#880000")
		brks=seq(0, max(sortedScores, na.rm=TRUE), length.out=30)
	}	
	
	png(paste("AVG_", n, "_heatmap.png", sep=""), width=4000, height=3000, res=300)
	aheatmap(sortedScores,
					 color=cols,
					 breaks=brks,
					 scale="none",
					 Rowv=NA, Colv=NA, 
					 cexRow=0,
					 cexCol=0)
	dev.off()
	
#	system(paste("composite ", paste("AVG_", n, "_heatmap.png", sep=""), " -compose Multiply chrom_contours_thick.png ", paste("AVG_", n, "_heatmap_contours_thick.png", sep=""), sep=""))
	system(paste("composite ", paste("AVG_", n, "_heatmap.png", sep=""), " -compose Multiply chrom_contours.png ", paste("AVG_", n, "_heatmap_contours.png", sep=""), sep=""))
}

# For the individual cancer types
for (n in names(results)) {
	for (sc in c("combined.score", "ts.score", "og.score")) {
		sortedScores = matrix(c(results[[n]]$prioritize.combined[sortedGenes, sc], rep(0, times=(139^2)-length(sortedGenes))), nrow=139, byrow=FALSE)[, 1:138]
		if (sc == "combined.score") {
			cols = colorpanel(29, low="#034b87", mid="white", high="#880000")
			ext = max(abs(min(sortedScores, na.rm=TRUE)), max(sortedScores, na.rm=TRUE))
			#brks=seq(-1*ext, ext, length.out=30)
			brks = 0
		} else if (sc == "ts.score") {
			cols = colorpanel(29, low="white", high="#034b87")
			brks=seq(0, max(sortedScores, na.rm=TRUE), length.out=30)
		} else {
			cols = colorpanel(29, low="white", high="#880000")
			brks=seq(0, max(sortedScores, na.rm=TRUE), length.out=30)
		}	
	
		png(paste(n, "_", sc, "_heatmap.png", sep=""), width=4000, height=3000, res=300)
		aheatmap(sortedScores,
						 color=cols,
					 	 breaks=brks,
					 	 scale="none",
					 	 Rowv=NA, Colv=NA, 
					 	 cexRow=0, cexCol=0)
		dev.off()
		
		system(paste("composite ", paste(n, "_", sc, "_heatmap.png", sep=""), " -compose Multiply chrom_contours.png ", paste(n, "_", sc, "_heatmap.png", sep=""), sep=""))
	}
}

