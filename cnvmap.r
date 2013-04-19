# Queries Biomart (NCBI36) to obtain the chromosomal location for each input gene.
# genes has to be a character vector with HGNC symbols
# Returns a data.frame with the chromosome, start and end positions
getGeneLocs = function(genes, warningsfile="warnings.txt", genome="hg18") {
  require(biomaRt)
  
  # List all Biomarts available through the corresponding archive
  # listMarts(host="may2009.archive.ensembl.org",path="/biomart/martservice")
  
  # Use the latest Biomart that is based on NCBI36/hg18
  mart = useMart(host="may2009.archive.ensembl.org", path="/biomart/martservice", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  if (genome == "hg19") {
  	mart = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  }
  
  # Map symbols that cannot be found in Ensembl to synonyms
  #hgnc.mart = useMart(host="www.genenames.org", path="/biomart/martservice", biomart="g4public", dataset="hgnc")
  # Find out the gene locations
  res = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters=c("hgnc_symbol"), values=genes, mart=mart)
  res = subset(res, chromosome_name %in% c(1:22, "X", "Y"))
  
  wf = NULL
  
  # Gene symbols for which no location could be found
  #notFound = setdiff(genes, res[, "hgnc_symbol"])
  #for (g in notFound) {
  #  # Map the current symbol to its synonyms
  #  syn = getBM(attributes=c("gd_prev_sym"), filters=c("gd_app_sym"), values=g, mart=hgnc.mart)
  #  if (nrow(syn) == 1) {
  #    # Get a vector with synonyms
  #    syn = unlist(strsplit(syn[1, 1], ", "))
  #    for (gs in syn) {
  #    	# Check whether Ensembl knows this symbol
  #      loc = getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters=c("hgnc_symbol"), values=gs, mart=mart)
  #      if (nrow(loc) > 0) {
  #        if (is.null(wf)) {
  #          wf = file(warningsfile, "a")
  #        }
  #        cat(paste("Mapped symbol ", g, " to ", loc[, "hgnc_symbol"], ", ", loc[, "chromosome_name"], ", ", loc[, "start_position"], ", ", loc[, "end_position"], "!\n", sep=""), file=wf)
  #        # Use the original symbol and add the location to the result
  #        loc[1, "hgnc_symbol"] = g
  #        res = rbind(res, loc)
  #        break
  #      }
  #    }
  #  } else {
  #    if (is.null(wf)) {
  #      wf = file(warningsfile, "a")
  #    }
  #    cat(paste("Could not map ", g, " to chromosomal location!\n", sep=""), file=wf)
  #  }
  #}
  
  if (length(unique(res$hgnc_symbol)) < length(res$hgnc_symbol)) {
  	if (is.null(wf)) {
  		wf = file(warningsfile, "a")
    }
    cat(paste("Found several genomic locations for genes: ", res$hgnc_symbol[which(duplicated(res$hgnc_symbol))], "\n", sep=""), file=wf)
    cat("Using the first locations!", file=wf)
    res = res[!duplicated(res$hgnc_symbol), ]
  }
  
  if (!is.null(wf)) {
    close(wf)
  }
  
  rownames(res) = res[, "hgnc_symbol"]
  if (length(which(res[, "chromosome_name"] == "X")) > 0) {
    res[which(res[, "chromosome_name"] == "X"), "chromosome_name"] = 23
  }
  if (length(which(res[, "chromosome_name"] == "Y")) > 0) {
    res[which(res[, "chromosome_name"] == "Y"), "chromosome_name"] = 24
  }
  return(res)
}

# Find the correct copy number segment and the associated copy number.
# min is the start of the location of interest
# min is the end of the location of interest
# segments is a matrix of start (column 1) and end (column 2) positions, and the
# copy number value (column 3) of the segments; all segments should be on the
# same chromosome as the position of interest
# If a gene spans several segments, the mean copy number of all these segments
# is returned.
getCopyNumberValue = function(start, stop, segments) {
  if (!is.matrix(segments)) {
    # Segments is a vector with only one segment, so add one artificial to make 
    # it a matrix. 
    segments = rbind(c(segments[1], 0, 0, 0), segments)
  }

  # Find the start segment
  # Look for the biggest segment start that is still smaller than the gene start
  startSeg = 0
  startSeg = max(which(segments[, 2] <= start))
  # If the end of that segment is larger than the start, the gene doesn't touch
  # that segment.
  startWithin = TRUE
  if (is.infinite(startSeg) || segments[startSeg, 3] < start)
    startWithin = FALSE
 
  # and the end segment
  # Look for the smallest segment end that is still bigger than the gene end
  endSeg = 0
  endSeg = min(which(segments[, 3] >= stop))
  # If the start of that segment is larger than the end, the gene doesn't touch
  # that segment.
  endWithin = TRUE
  if (is.infinite(endSeg) || segments[endSeg, 2] > stop)
    endWithin = FALSE

 
  if (is.infinite(startSeg) && !is.infinite(endSeg) || (startSeg > endSeg)) {
    # If no start segment but an end segment was found, the gene begins before the first one.
    # Use the first segment in this case
    startSeg = 1
  } else if (!is.infinite(startSeg) && is.infinite(endSeg)) {
    # If no end segment but a start segment was found, the gene ends behind the last one.
    # Use the last segment in this case
    endSeg = nrow(segments)
  }

  case = 0
  if (is.infinite(startSeg) && is.infinite(endSeg)) {
    # The gene doesn't touch any segments
    cnvalue = NA
    case = 1
  } else if (startSeg == endSeg) {
    # The gene lies withing one segment
    cnvalue = segments[startSeg, 4]
    case = 2
  } else if (!startWithin) {
    # The gene's start lies between segments
    # -> ignore the start segment and average the remaining
    cnvalue = mean(segments[(startSeg+1):endSeg, 4])
    case = 3
  } else if (!endWithin) {
    # The gene's end lies between segments
    # -> ignore the end segment and average the remaining
    cnvalue = mean(segments[startSeg:(endSeg-1), 4])
    case = 4
  } else {
    # The gene spans different segments
    # Take the average of all copy number values of these segments
    cnvalue = mean(segments[startSeg:endSeg, 4])
    case = 5
  }

  return(list(cnvalue=as.double(cnvalue), case=case))
}

# Maps all input genes to the copy number value of the corresponding segment.
# genes has to be a character vector containing HGNC symbols
# cnData should be a matrix with segmented copy number variation data; formatted
# as CBS output
# cases is a boolean value; if set to TRUE, an additional matrix containing 
# information on how each CN value was determined is returned 
mapGenes2CN = function (genes, cnData, cases=FALSE, genome="hg18") {
  # Map all genes to their chromosomal location using biomaRt
  genes2loc = getGeneLocs(genes, genome=genome)
  # List of samples
  allSamples = unique(cnData[, 1])
  
  # Matrix that will contain the copy number for each gene in each sample
  cn.matrix = matrix(NA, nrow=length(genes), ncol=length(allSamples))
  rownames(cn.matrix) = genes
  colnames(cn.matrix) = allSamples
  
  # If the case matrix is to be reported
  if (cases) {
    cases = matrix(-1, nrow=length(genes), ncol=length(allSamples))
    rownames(cases) = genes
    colnames(cases) = allSamples
  }
  
  # Go through all samples
  for (samp in allSamples) {
    # Get all segments with their copy number for the current sample
    samp.segs = as.matrix(cnData[which(cnData[, 1] == samp), c(2, 3, 4, 6)])
    # Go through the chromosomes
    for (chrom in unique(genes2loc[, "chromosome_name"])) {
      # Get the segments for that particular chromosome
      segs = samp.segs[which(samp.segs[, 1] == chrom), ]
      # Go through all genes on that chromosome
      for (gene in rownames(genes2loc[which(genes2loc[, "chromosome_name"] == chrom), ])) {
      	# Get the copy numberk
        res = getCopyNumberValue(as.integer(genes2loc[gene, "start_position"]), as.integer(genes2loc[gene, "end_position"]), segs)
        cn.matrix[gene, samp] = res[["cnvalue"]]
        # Remember the case if that information is wanted
        if (cases) {
          cases[gene, samp] = res[["case"]]
        }
      }
    }
  }
  
  res = list(cn.values=cn.matrix)
  if (cases) {
    res[["cn.cases"]] = cases
  }
  return(res)
}

# Compares aCGH data of tumor and normal samples.
# The function will run paired Wilcoxon tests. Pairing of samples
# is determined by matching sample names in the tumors and normals matrices.
# The function uses doWilcox() as defined for the analysis of methylation data!
# tumors is an aCGH data matrix with samples in the columns and genes in the rows
# normals is a methylation data matrix with samples in the columns and genes in the rows
# Returns a named list with results from paired tests.
runCGHComp = function(tumors, normals) {
  # Samples from tumors that have a matched normal
  tumors.matchedsamples = intersect(colnames(tumors), colnames(normals))
  # Rename normal samples to reflect the state
  colnames(normals) = paste(colnames(normals), "normal", sep="_")

  tumors.groups = c(rep(1, times=ncol(tumors)), rep(2, times=ncol(normals)))
  names(tumors.groups) = c(colnames(tumors), colnames(normals))

  # Run paired Wilcoxon tests
  if (nrow(tumors) > 1) {
    inpMat = cbind(tumors[, tumors.matchedsamples], normals)
  } else {
    inpMat = matrix(c(tumors[, tumors.matchedsamples], normals), nrow=1, byrow=TRUE)
    rownames(inpMat) = rownames(tumors)
    colnames(inpMat) = c(tumors.matchedsamples, colnames(normals))
  }
  pairedwilcox.p = doWilcox(inpMat=inpMat, matchedSamples=tumors.matchedsamples)
 
  return(list(paired=pairedwilcox.p))
}

# Calculates correlation between copy number and expression of the corresponding genes.
# Correlation is calculated using all samples contained in both data matrices
# cnv.data should be a matrix with genes in rows and samples in columns
# exprs.data should be a matrix with probes in rows and samples in columns
# genes should be a character vector giving the genes to test
corCnvExprs = function(cnv.data, exprs.data, genes) {
  # All genes contained in the expression data
  exprs.genes = rownames(exprs.data)
  commonSamples = intersect(colnames(cnv.data), colnames(exprs.data))
	
  cors = c()
  n = c()
  for (x in intersect(genes, exprs.genes)) {
    cors = c(cors, cor(exprs.data[x, commonSamples], cnv.data[x, commonSamples], method="spearman", use="pairwise.complete.obs"))
    n = c(n, x)
  }
  names(cors) = n
  
  return(cors)
}

# For each gene, summarizes evidence for significant copy number alterations and the correlation with gene expression. 
# cnv.wilcox is a vector with p-values (preferably multiple-testing corrected) for copy number difference between tumor and normal
# cnv.diff is a vector giving the value mean(cnv.tumor) - mean(cnv.normal) for each gene
# cors is a vector containing the correlation between copy number and expression for each gene
# genes is a character vector of genes to consider
# wilcox.cutoff is the significance cut-off, probes with cnv.wilcox > wilcox.cutoff will be considered significant
# diff.cutoff genes need to have a copy number difference smaller or greater than this cut-off to be considered significant
# regulation either "down" or "up" for finding genes that are lost or gained, respectively
# cor.min genes need to have a correlation larger or equal than this cut-off to be considered significant
# cor.max genes need to have a correlation smaller or equal than this cut-off to be considered significant
summarizeCnvExprs = function(cnv.wilcox, cnv.diff, cors, genes, wilcox.cutoff=0.05, diff.cutoff=0.1, regulation=c("down", "up"), cor.min=0.1, cor.max=1.0) { 
  regulation = match.arg(regulation)
  
  # Get the correct comparison function
  # If we want to find genes with higher copy number in tumors get the greaterThan function
  # If we want to find genes with lower copy number in tumors, get the smallerThan function
  compare = switch(regulation, down=smallerThan, up=greaterThan)
  
  gene.scores = rep(0.0, length(genes))
  names(gene.scores) = genes
  
  for (gene in genes) {
    # Get p-value for the copy number difference for that gene
    significant = cnv.wilcox[gene] < wilcox.cutoff
    # Check whether the copy number is correlated with expression
    correlated = (cors[gene] >= cor.min && cors[gene] <= cor.max)
    # Get the actual methylation differences
    difference = compare(cnv.diff[gene], diff.cutoff)
    # Does the gene pass all filters?
    gene.scores[gene] = as.integer(significant && correlated && difference) 
  }
  
  return(gene.scores)
}


doCnvAnalysis = function(tumors, normals, genes, exprs, samples=NULL, wilcox.cutoff=0.05, diff.cutoff=-0.1, regulation=c("down", "up"), cor.min=0.1, cor.max=1.0) {
  regulation = match.arg(regulation)
  
  # Get the correct comparison function
  # If we want to find genes with higher copy number in tumors get the greaterThan function
  # If we want to find genes with lower copy number in tumors, get the smallerThan function
  compare = switch(regulation, down=smallerThan, up=greaterThan)
  
  if (is.null(samples)) {
    tumor.samples = colnames(tumors)
    normal.samples = colnames(normals)
  } else {
    tumor.samples = intersect(samples, colnames(tumors))
    normal.samples = intersect(samples, colnames(normals))
  }
  
  selected.genes = intersect(genes, intersect(rownames(tumors), rownames(normals)))
  
  # Calculate the difference in copy number between tumors and normals
  # tumor - normal
  # --> diff < 0 implies that the gene has higher mean copy number values in normals than in tumors
  # --> diff > 0 implies that the gene has lower mean copy number values in normals than in tumors
  cgh.mean.tumor = apply(tumors[selected.genes, tumor.samples], 1, mean, na.rm=TRUE)
  cgh.mean.normal = apply(normals[selected.genes, normal.samples], 1, mean, na.rm=TRUE)
  cgh.mean.diff = cgh.mean.tumor - cgh.mean.normal
  selected.genes = names(cgh.mean.diff[which(compare(cgh.mean.diff, diff.cutoff))])
  
  cgh.cors = list()
  if (length(selected.genes) > 0) {
    temp = matrix(tumors[selected.genes, tumor.samples], nrow=length(selected.genes), byrow=TRUE)
    rownames(temp) = selected.genes
    colnames(temp) = tumor.samples
    cgh.cors = list(paired=corCnvExprs(temp, exprs, selected.genes))
    selected.genes = names(cgh.cors$paired[which(cgh.cors$paired >= cor.min & cgh.cors$paired <= cor.max)])
  }
  
  cgh.wilcox.bh = list()
  cgh.wilcox = list()
  if (length(selected.genes) > 0) {
    temp = matrix(tumors[selected.genes, tumor.samples], nrow=length(selected.genes), byrow=TRUE)
    rownames(temp) = selected.genes
    colnames(temp) = tumor.samples
    temp.n = matrix(normals[selected.genes, normal.samples], nrow=length(selected.genes), byrow=TRUE)
    rownames(temp.n) = selected.genes
    colnames(temp.n) = normal.samples
    cgh.wilcox = runCGHComp(temp, temp.n)
    # Multiple testing correction
    cgh.wilcox.bh = list(paired=p.adjust(cgh.wilcox$paired, method="BH"))
    selected.genes = names(cgh.wilcox.bh$paired[which(cgh.wilcox.bh$paired < wilcox.cutoff)])
  }
  
  gene.scores = rep(0, length(genes))
  names(gene.scores) = genes
  gene.scores[selected.genes] = 1
                                          
  return(list(scores=gene.scores, cors=cgh.cors, wilcox=cgh.wilcox, corrected=cgh.wilcox.bh, diffs=cgh.mean.diff))
}

