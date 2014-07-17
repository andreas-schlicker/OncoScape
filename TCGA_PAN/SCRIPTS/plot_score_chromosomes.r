# Load libraries 
library(biomaRt)
library(Gviz)
library(GenomicRanges)

# Prioritization results
load("~/bulk/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata")

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

# For each score, plot each chromosome across all cancer types
sortedCancers = sort(names(results))
for (sc in c("combined.score", "ts.score", "og.score", "ts.affected.rel", "og.affected.rel", "og.ts.affected.rel")) {
	for (chrom in 1:24) {
		genes = subset(geneLoc, chromosome_name==chrom)
		
		tsScore = length(grep("ts.", sc)) > 0 && sc != "og.ts.aggected.rel"
		
		if (sc != "og.ts.affected.rel") {
			scoresAvg = results.combined[genes[, "hgnc_symbol"], sc] * ifelse(tsScore, -1, 1)
		} else {
			scoresAvg = results.combined[genes[, "hgnc_symbol"], "og.affected.rel"] - results.combined[genes[, "hgnc_symbol"], "ts.affected.rel"]
		}
		
		tracks = c(list(AVG=DataTrack(GRanges(seqnames=as.character(chrom),
													 								ranges=IRanges(start=genes[, "start_position"],
													 												 			 end=genes[, "end_position"]),
													 								scoresAvg,
													 								strand="*"),
													 				name="AVG", background.title="gray85", background.panel="gray85")),
													 				
							 lapply(1:length(sortedCancers), function(x) {
							 	 if (sc != "og.ts.affected.rel") {
							 	 	 scores = results[[sortedCancers[x]]]$prioritize.combined[genes[, "hgnc_symbol"], sc] * ifelse(tsScore, -1, 1)
							 	 } else {
							 	 	 scores = results[[sortedCancers[x]]]$prioritize.combined[genes[, "hgnc_symbol"], "og.affected.rel"] -
							 	 	 					results[[sortedCancers[x]]]$prioritize.combined[genes[, "hgnc_symbol"], "ts.affected.rel"]
							 	 }
							 	 
							 	 DataTrack(GRanges(seqnames=as.character(chrom),
													 				 ranges=IRanges(start=genes[, "start_position"],
													 									 			end=genes[, "end_position"]),
													 				 scores,
							 										 strand="*"),
							 							name=sortedCancers[x],
							 							background.title=ifelse(x %% 2 == 1, "white", "gray85"),
							 							background.panel=ifelse(x %% 2 == 1, "white", "gray85"))}))
							 							
		png(paste("chromosome", chrom, "_", sc, ".png", sep=""), width=4000, height=3000, res=300)							 							
		plotTracks(tracks, type=c("h", "g"), ylim=c(min(unlist(lapply(tracks, function(x) { min(x@data) }))), 
																								max(unlist(lapply(tracks, function(x) { max(x@data) })))), 
							 col.title="black", col.axis="black", cex.axis=0.75, fontface="bold", 
							 main=paste("Chromosome ", chrom, sep=""))
		dev.off()

		
		if (sc == "og.ts.affected.rel") {
			limits = c(-1, 1)
		} else {
			limits = range(results.combined[, sc] * ifelse(tsScore, -1, 1))
		}

		# Plot the average scores for each chromosome and score type		
		png(paste("chromosome", chrom, "_all_avg_", sc, ".png", sep=""), width=3000, height=1500, res=300)
		genes = subset(geneLoc, chromosome_name==chrom)
		plotTracks(DataTrack(GRanges(seqnames=as.character(chrom),
													 		 	 ranges=IRanges(start=genes[, "start_position"],
													 											end=genes[, "end_position"]),
													 			 scoresAvg,
													 			 strand="*"),
												 #name=paste("Chr", chrom, sep=""), 
												 background.title=ifelse(chrom %% 2 == 1, "white", "gray85"),
												 background.panel=ifelse(chrom %% 2 == 1, "white", "gray85")),
						   type=c("h", "g"), ylim=limits, 
							 col.title="black", col.axis="black", cex.axis=0.75, fontface="bold",
							 panel.only=TRUE)
		dev.off()
	}
}			


# Combine them into one plot 
#montage chromosome1_all_avg_combined.score.png  chromosome2_all_avg_combined.score.png  chromosome3_all_avg_combined.score.png  chromosome4_all_avg_combined.score.png  chromosome5_all_avg_combined.score.png  chromosome6_all_avg_combined.score.png  chromosome7_all_avg_combined.score.png  chromosome8_all_avg_combined.score.png  chromosome9_all_avg_combined.score.png  chromosome10_all_avg_combined.score.png  chromosome11_all_avg_combined.score.png  chromosome12_all_avg_combined.score.png  chromosome13_all_avg_combined.score.png  chromosome14_all_avg_combined.score.png  chromosome15_all_avg_combined.score.png  chromosome16_all_avg_combined.score.png  chromosome17_all_avg_combined.score.png  chromosome18_all_avg_combined.score.png  chromosome19_all_avg_combined.score.png  chromosome20_all_avg_combined.score.png  chromosome21_all_avg_combined.score.png  chromosome22_all_avg_combined.score.png  chromosome23_all_avg_combined.score.png  chromosome24_all_avg_combined.score.png  -background white -tile 5x5 -geometry +0+0 PNG32:chromosome_all_avg_combined.score.png
#montage chromosome1_all_avg_ts.score.png  chromosome2_all_avg_ts.score.png  chromosome3_all_avg_ts.score.png  chromosome4_all_avg_ts.score.png  chromosome5_all_avg_ts.score.png  chromosome6_all_avg_ts.score.png  chromosome7_all_avg_ts.score.png  chromosome8_all_avg_ts.score.png  chromosome9_all_avg_ts.score.png  chromosome10_all_avg_ts.score.png  chromosome11_all_avg_ts.score.png  chromosome12_all_avg_ts.score.png  chromosome13_all_avg_ts.score.png  chromosome14_all_avg_ts.score.png  chromosome15_all_avg_ts.score.png  chromosome16_all_avg_ts.score.png  chromosome17_all_avg_ts.score.png  chromosome18_all_avg_ts.score.png  chromosome19_all_avg_ts.score.png  chromosome20_all_avg_ts.score.png  chromosome21_all_avg_ts.score.png  chromosome22_all_avg_ts.score.png  chromosome23_all_avg_ts.score.png  chromosome24_all_avg_ts.score.png  -background white -tile 5x5 -geometry +0+0 PNG32:chromosome_all_avg_ts.score.png
#montage chromosome1_all_avg_og.score.png  chromosome2_all_avg_og.score.png  chromosome3_all_avg_og.score.png  chromosome4_all_avg_og.score.png  chromosome5_all_avg_og.score.png  chromosome6_all_avg_og.score.png  chromosome7_all_avg_og.score.png  chromosome8_all_avg_og.score.png  chromosome9_all_avg_og.score.png  chromosome10_all_avg_og.score.png  chromosome11_all_avg_og.score.png  chromosome12_all_avg_og.score.png  chromosome13_all_avg_og.score.png  chromosome14_all_avg_og.score.png  chromosome15_all_avg_og.score.png  chromosome16_all_avg_og.score.png  chromosome17_all_avg_og.score.png  chromosome18_all_avg_og.score.png  chromosome19_all_avg_og.score.png  chromosome20_all_avg_og.score.png  chromosome21_all_avg_og.score.png  chromosome22_all_avg_og.score.png  chromosome23_all_avg_og.score.png  chromosome24_all_avg_og.score.png  -background white -tile 5x5 -geometry +0+0 PNG32:chromosome_all_avg_og.score.png
#montage chromosome1_all_avg_ts.affected.rel.png  chromosome2_all_avg_ts.affected.rel.png  chromosome3_all_avg_ts.affected.rel.png  chromosome4_all_avg_ts.affected.rel.png  chromosome5_all_avg_ts.affected.rel.png  chromosome6_all_avg_ts.affected.rel.png  chromosome7_all_avg_ts.affected.rel.png  chromosome8_all_avg_ts.affected.rel.png  chromosome9_all_avg_ts.affected.rel.png  chromosome10_all_avg_ts.affected.rel.png  chromosome11_all_avg_ts.affected.rel.png  chromosome12_all_avg_ts.affected.rel.png  chromosome13_all_avg_ts.affected.rel.png  chromosome14_all_avg_ts.affected.rel.png  chromosome15_all_avg_ts.affected.rel.png  chromosome16_all_avg_ts.affected.rel.png  chromosome17_all_avg_ts.affected.rel.png  chromosome18_all_avg_ts.affected.rel.png  chromosome19_all_avg_ts.affected.rel.png  chromosome20_all_avg_ts.affected.rel.png  chromosome21_all_avg_ts.affected.rel.png  chromosome22_all_avg_ts.affected.rel.png  chromosome23_all_avg_ts.affected.rel.png  chromosome24_all_avg_ts.affected.rel.png  -background white -tile 5x5 -geometry +0+0 PNG32:chromosome_all_avg_ts.affected.rel.png
#montage chromosome1_all_avg_og.affected.rel.png  chromosome2_all_avg_og.affected.rel.png  chromosome3_all_avg_og.affected.rel.png  chromosome4_all_avg_og.affected.rel.png  chromosome5_all_avg_og.affected.rel.png  chromosome6_all_avg_og.affected.rel.png  chromosome7_all_avg_og.affected.rel.png  chromosome8_all_avg_og.affected.rel.png  chromosome9_all_avg_og.affected.rel.png  chromosome10_all_avg_og.affected.rel.png  chromosome11_all_avg_og.affected.rel.png  chromosome12_all_avg_og.affected.rel.png  chromosome13_all_avg_og.affected.rel.png  chromosome14_all_avg_og.affected.rel.png  chromosome15_all_avg_og.affected.rel.png  chromosome16_all_avg_og.affected.rel.png  chromosome17_all_avg_og.affected.rel.png  chromosome18_all_avg_og.affected.rel.png  chromosome19_all_avg_og.affected.rel.png  chromosome20_all_avg_og.affected.rel.png  chromosome21_all_avg_og.affected.rel.png  chromosome22_all_avg_og.affected.rel.png  chromosome23_all_avg_og.affected.rel.png  chromosome24_all_avg_og.affected.rel.png  -background white -tile 5x5 -geometry +0+0 PNG32:chromosome_all_avg_og.affected.rel.png
#montage chromosome1_all_avg_og.ts.affected.rel.png  chromosome2_all_avg_og.ts.affected.rel.png  chromosome3_all_avg_og.ts.affected.rel.png  chromosome4_all_avg_og.ts.affected.rel.png  chromosome5_all_avg_og.ts.affected.rel.png  chromosome6_all_avg_og.ts.affected.rel.png  chromosome7_all_avg_og.ts.affected.rel.png  chromosome8_all_avg_og.ts.affected.rel.png  chromosome9_all_avg_og.ts.affected.rel.png  chromosome10_all_avg_og.ts.affected.rel.png  chromosome11_all_avg_og.ts.affected.rel.png  chromosome12_all_avg_og.ts.affected.rel.png  chromosome13_all_avg_og.ts.affected.rel.png  chromosome14_all_avg_og.ts.affected.rel.png  chromosome15_all_avg_og.ts.affected.rel.png  chromosome16_all_avg_og.ts.affected.rel.png  chromosome17_all_avg_og.ts.affected.rel.png  chromosome18_all_avg_og.ts.affected.rel.png  chromosome19_all_avg_og.ts.affected.rel.png  chromosome20_all_avg_og.ts.affected.rel.png  chromosome21_all_avg_og.ts.affected.rel.png  chromosome22_all_avg_og.ts.affected.rel.png  chromosome23_all_avg_og.ts.affected.rel.png  chromosome24_all_avg_og.ts.affected.rel.png  -background white -tile 5x5 -geometry +0+0 PNG32:chromosome_all_avg_og.ts.affected.rel.png


## Combined different scores per chromosome and cancer type
scores = c("combined.score", "og.score", "ts.score")
names(scores) = c("Combined score", "OG score", "TS score")

for (canc in sortedCancers) {
	for (chrom in 1:24) {
		genes = subset(geneLoc, chromosome_name==chrom)
		tracks = c(lapply(1:length(scores), function(x) { 
											DataTrack(GRanges(seqnames=as.character(chrom),
															 				  ranges=IRanges(start=genes[, "start_position"],
															 												 end=genes[, "end_position"]),
															 					results[[canc]]$prioritize.combined[genes[, "hgnc_symbol"], scores[x]] * 
															 						ifelse(length(grep("ts.", scores[x])) > 0, -1, 1),
															 					strand="*"),
															 	name=names(scores)[x], 
															 	background.title=ifelse(x %% 2 == 1, "white", "gray85"),
									 							background.panel=ifelse(x %% 2 == 1, "white", "gray85"))}),
		
								DataTrack(GRanges(seqnames=as.character(chrom),
															 		ranges=IRanges(start=rep(genes[, "start_position"], times=2),
															 									 end=rep(genes[, "end_position"], times=2)),
															 		c(results[[canc]]$prioritize.combined[genes[, "hgnc_symbol"], "og.affected.rel"], 
															 			-1*results[[canc]]$prioritize.combined[genes[, "hgnc_symbol"], "ts.affected.rel"]),
									 								strand="*"),
									 				name="% affected samples",
									 				background.title="gray85",
									 				background.panel="gray85"))
		
		png(paste(canc, "_chromosome", chrom, ".png", sep=""), width=4000, height=3000, res=300)
		plotTracks(tracks, type=c("h", "g"),
							 ylim=c(min(unlist(lapply(tracks, function(x) { min(x@data) }))), 
											max(unlist(lapply(tracks, function(x) { max(x@data) })))),
							 col.title="black", col.axis="black", cex.axis=0.75, fontface="bold", 
							 main=paste("Chromosome ", chrom, sep=""))
		dev.off()
	}
}

