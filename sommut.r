##' Creates the mutation ID from the given indexes.
##' If "reference" is set, "tumor1" and "tumor2" have to be set as well.
##' Only the allele that is different from the reference is added to the ID. 
##' @param data a mutation row
##' @param indexes vector with indexes to use
##' @param reference column index for the reference allele; default: NULL
##' @param tumor1 column index for the first tumor allele; default: NULL
##' @param tumor2 column index for the second tumor allele; default: NULL
##' @return mutation ID
##' @author Andreas Schlicker
mutationId = function(data, indexes, reference=NULL, tumor1=NULL, tumor2=NULL) {
	if (!is.null(reference)) {
		if (data[reference] == data[tumor1]) {
			indexes = c(indexes, tumor2)
		} else {
			indexes = c(indexes, tumor1)
		}
	}
	
	paste(data[indexes], collapse="_")
}

##' Analyzes the mutation data.
##' The method calculates oncogene and tumor suppressor gene scores based on the
##' 20/20 rule by Vogelstein et al. (doi: 10.1126/science.1235122). 
##' @param mutations mutation matrix using TCGA format
##' @param genes character vector with gene IDs; default: NULL (use all genes)
##' @param ignore character vector with mutation types to ignore; default: "silent"
##' @param samples vector of samples to restrict the analysis to; default: NULL
##' @param filtercol index of the column that should be used for filtering; default: NULL
##' (no filtering)
##' @param filter threshold for filtering mutations; only mutations with score >= "filter"
##' will be taken into account; default: 1
##' @param genecol index of the gene id column; default: 1
##' @param typecol index of the mutation type column; default: 9
##' @param samplecol index of the sample id column; default: 16
##' @param chromcol index of the chromosome column; default: 5
##' @param startcol index of the column with start position; default: 6
##' @param endcol index of the column with end position; default: 7
##' @param reference index of the column with reference allele; default: 11
##' @param tumor1 index of the column with tumor allele 1; default: 12
##' @param tumor2 index of the column with tumor allele 2; default: 13
##' @return matrix with genes in rows and values in columns
##' @author Andreas Schlicker
doMutationAnalysis = function(mutations, genes=NULL,  
							  ignore=c("Silent"), samples=NULL,
							  filtercol=NULL, filter=1,
                       		  genecol=1, typecol=9, samplecol=16, 
                       		  chromcol=5, startcol=6, endcol=7,
							  reference=11, tumor1=12, tumor2=13) {
				   
	if (is.null(genes)) { 
		genes = unique(mutations[, genecol])
	}
	
	# Filter samples
	if (!is.null(samples)) {
		mutations = mutations[which(mutations[, samplecol] %in% samples), ]
	}
	
	if (!is.null(ignore) && length(ignore) > 0) {
		mutations = mutations[which(!(mutations[, typecol] %in% ignore)), ]
	}
	
	if (!is.null(filtercol)) {
		mutations = mutations[which(mutations[, filtercol] >= filter), ]
	}
	
	res = matrix(NA, nrow=length(genes), ncol=8)
	rownames(res) = genes
	colnames(res) = c("og", "ts", "total.mutations", "unique.mutations", "mutated.samples", "total.samples", "ts.mutations", "og.mutations")
	res[, "total.samples"] = length(unique(mutations[, samplecol]))
	
	for (gene in genes) {
		# Get all rows with mutations in the current gene that shouldn't be ignored
		muts = mutations[which(mutations[, genecol] == gene), ]
		
		missense = table(apply(muts[which(muts[, typecol] == "Missense_Mutation"), ], 1, 
						 function(x) { mutationId(x, indexes=c(chromcol, startcol, endcol)) } ))
		inframe_del = table(apply(muts[which(muts[, typecol] == "In_Frame_Del"), ], 1, 
						function(x) { mutationId(x, indexes=c(chromcol, startcol, endcol)) } ))
		inframe_ins = table(apply(muts[which(muts[, typecol] == "In_Frame_Ins"), ], 1, 
						function(x) { mutationId(x, indexes=c(chromcol, startcol, endcol), reference, tumor1, tumor2) } ))
		
		res[gene, "og"] = 1 - ((length(missense) + length(inframe_del) + length(inframe_ins)) / nrow(muts))
		#res[gene, "og"] = 1 - ((length(missense) + length(inframe_del) + length(inframe_ins)) / length(which(muts[, typecol] %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"))))
		
		ts.muts = c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site")
		res[gene, "ts"] = nrow(muts[which(muts[, typecol] %in% ts.muts), ]) / nrow(muts)
		
		# Count total number of mutations
		res[gene, "total.mutations"] = nrow(muts)
		# Count number of distinct mutations
		# Distinct is defined as unique combination of chromosome, start and end positions
		res[gene, "unique.mutations"] = length(unique(apply(muts, 1, function(x) { mutationId(x, c(chromcol, startcol, endcol)) })))
		# Count the number of distinct samples with mutations in that gene
		res[gene, "mutated.samples"] = length(unique(muts[, samplecol]))
		res[gene, "og.mutations"] = sum(missense) + sum(inframe_del) + sum(inframe_ins)
		res[gene, "ts.mutations"] = res[gene, "ts"] * nrow(muts)		
	}
	
	res
}

##' Summarizes mutation scores. Genes with an og.score > og.cutoff && ts.score < ts.cutoff.low
##' are classified as oncogenes. Genes with a ts.score > ts.cutoff || og.score > og.cutoff && ts.score > ts.cutoff.low
##' are classified as tumor suppressors.
##' @param mutations mutation matrix using TCGA format
##' @param mut.analysis output of doMutationAnalysis
##' @param genes vector of gene IDs to test; default: NULL (test all genes)
##' @param samples vector of sample IDs to consider; default: NULL (all samples)
##' @param og.cutoff cut-off for the oncogene score; default: 0.2
##' @param ts.cutoff cut-off for the tumor suppressor score; default: 0.2
##' @param ts.cutoff.log cut-off for the low tumor suppressor score; default: 0.05
##' @param score either "ts" or "og" to score tumor suppressors or oncogenes, respectively
##' @param genecol index of the gene id column; default: 1
##' @param typecol index of the mutation type column; default: 9
##' @param samplecol index of the sample id column; default: 16
##' @return named list with scores for genes ("scores"), number of affected samples ("summary") and the lists of affected samples ("samples")
##' @author Andreas Schlicker
summarizeMutations = function(mutations, mut.analysis, 
							  genes=NULL, samples=NULL,
							  og.cutoff=0.2, ts.cutoff=0.2, ts.cutoff.low=0.05, 
							  score=c("ts", "og"),
							  genecol=1, typecol=9, samplecol=16) {
	score = match.arg(score)
	
  	if (is.null(genes)) {
		genes = rownames(mut.analysis)
  	} else {
		mutations = mutations[which(mutations[, genecol] %in% genes), ]
	}
	common = intersect(genes, mutations[, genecol])
	missing = setdiff(genes, mutations[, genecol])
	
	if (!is.null(samples)) {
		mutations = mutations[which(mutations[, samplecol] %in% samples), ]
	}
  
	if (score == "ts") {
		gene.scores = as.integer((mut.analysis[, "ts.mutations"] >= 5) & (mut.analysis[, "ts"] > ts.cutoff | (mut.analysis[, "og"] > og.cutoff & mut.analysis[, "ts"] > ts.cutoff.low))[common])
	} else {
		gene.scores = as.integer((mut.analysis[, "og.mutations"] >= 5) & (mut.analysis[, "og"] > og.cutoff & mut.analysis[, "ts"] < ts.cutoff.low)[common])
	}
	names(gene.scores) = common
	gene.scores[missing] = 0
	gene.scores = gene.scores[genes]
	
	affected.samples = mutationsAffectedSamples(genes, names(gene.scores)[gene.scores == 1], mutations,
												score, genecol, typecol, samplecol)
	list(scores=gene.scores, summary=affected.samples$summary, samples=affected.samples$samples)
}

##' Count the number of tumors affected by mutations in the given genes.
##' @param genes vector with gene IDs that will be contained in the results
##' @param test.genes vector with gene IDs that will be tested; default: genes
##' @param mutations the matrix with the mutation data
##' @param type either "ts" or "og" to get samples with tumor suppressor or oncogene mutations, respectively
##' @param genecol index of the gene id column; default: 1
##' @param typecol index of the mutation type column; default: 9
##' @param samplecol index of the sample id column; default: 16
##' @return named list with two components; "summary" is a matrix with absolute (1st column) 
##' and relative (2nd column) numbers of affected samples; "samples" is a named list with all 
##' samples affected by a change in this feature
##' @author Andreas Schlicker
mutationsAffectedSamples = function(genes, test.genes=genes,
									mutations, type=c("ts", "og"),
									genecol=1, typecol=9, samplecol=16) {  
	
	type = match.arg(type)
	# Filter mutations
	if (type == "ts") {
		mutations = mutations[which(mutations[, typecol] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site")), ]
	} else {
		mutations = mutations[which(mutations[, typecol] %in% c("Missense_Mutation", "In_Frame_Del", "In_Frame_Ins")), ]
	}	
	
	normFactor = length(unique(mutations[, samplecol]))
	
	# Which genes to test are in the mutation data?
	common = intersect(test.genes, unique(mutations[, genecol]))
	# Which genes should appear in the output but are missing?
	missing = setdiff(genes, common)
	
	samples = lapply(common, function(x) { unique(mutations[which(mutations[, genecol] == x), samplecol]) })
	names(samples) = common
	affected = unlist(lapply(samples, length))
	
	# Add all missing features and resort
	affected[missing] = 0
	affected = affected[genes]
	samples = samples[genes]
	names(samples) = genes
	
	list(summary=cbind(absolute=affected, relative=affected/normFactor), samples=samples)
}
