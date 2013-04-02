# Analyzes the mutation data.
# The method counts the number of mutations for each gene in the list. 
# Depending on the parameters, all mutations per gene are counted or the number
# of samples with mutations are counted for each gene.
# mutations the matrix with the mutation data
# genes a character vector containing the genes to be analyzed
# ignore character vector containing the type of mutations that should be 
# ignored, this vector is used as is; defaults to c("silent")
# genecol integer giving the index of the gene id column; defaults to 1
# typecol integer giving the index of the mutation type column; defaults to 9
# samplecol integer giving the index of the sample id column; defaults to 16
# chromcol integer giving the index of the chromosome column; defaults to 5
# startcol integer giving the index of the start column; defaults to 6
# stopcol integer giving the index of the stop column; defaults to 7
analyzeMuts = function(mutations, genes, ignore=c("Silent"), 
											 samples=NULL,
                       genecol=1, typecol=9, samplecol=16, 
                       chromcol=5, startcol=6, stopcol=7) {  
  res = matrix(0, nrow=length(genes), ncol=4)
  rownames(res) = genes
  colnames(res) = c("total mutations", "unique mutations", "mutated samples", "unique samples")
  
  locMutations = mutations
  if (!is.null(samples)) {
  	locMutations = mutations[which(mutations[, samplecol] %in% samples), ]
  }
  
  res[, 4] = length(unique(locMutations[, samplecol]))
  
  for (gene in genes) {
    # Get all rows with mutations in the current gene that shouldn't be ignored
    muts = locMutations[which(locMutations[, genecol] == gene & !(locMutations[, typecol] %in% ignore)), ]
    
    # Count total number of mutations
    res[gene, "total mutations"] = nrow(muts)
    # Count number of distinct mutations
    # Distinct is defined as unique combination of chromosome, start and end positions
    res[gene, "unique mutations"] = length(unique(apply(muts, 1, function(x) { paste(x[chromcol], x[startcol], x[stopcol], sep="_") })))
    # Count the number of distinct samples with mutations in that gene
    res[gene, "mutated samples"] = length(unique(muts[, samplecol]))
  }
  
  return(res)
}

# Summarizes mutation scores. 
# mutations matrix with genes in the rows and at least the following three columns: 
# columns 1 through 3 have to be total number of mutations, number of unique 
# mutations, and number of mutated samples, in this order; 
# e.g. the output of the "analyzeMuts" function
# genes a character vector containing the ids of the genes to analyze; ids have
# to map to row names in the "mutations" matrix
# ratio.mut.samples floating point number between 0 and 1 giving the minimum percentage of 
# samples that need to have a mutation in a certain gene; defaults to 0.01
# uniqueness floating point number between 0 and 1 giving the ratio of unique
# mutations required; the ratio is calculated as unique/total
# comparison either "down" or "up"; defines whether the ratio of unique mutations
# should be larger or smaller than "uniqueness"; oncogenes should have
# a high uniqueness value and tumor suppressor genes a low uniqueness value
# Returns a vector giving the score for each gene
summarizeMutations = function(mutations, genes, ratio.mut.samples=0.01, uniqueness=0.1, comparison=c("smaller", "larger")) {
  comparison = match.arg(comparison)
  
  # Get the correct comparison function
  # If we want to find genes with few unique mutations, get the smallerThan function
  # If we want to find genes with many unique mutations, get the greaterThan function
  compare = switch(comparison, smaller=smallerThan, larger=greaterThan)
  
  gene.scores = rep(0, length(genes))
  names(gene.scores) = genes
  
  for (gene in genes) {
    if (sum(is.na(mutations[gene, ])) == 0) {
      if (((mutations[gene, 3] / mutations[gene, 4]) > ratio.mut.samples) && compare(mutations[gene, 2] / mutations[gene, 1], uniqueness)) {
        gene.scores[gene] = 1
      }
    }
  }
  
  return(gene.scores)
}
