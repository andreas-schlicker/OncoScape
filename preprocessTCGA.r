# Renames the file by replacing the extract ID with the sample code
# path relative or absolute pathname
# pattern pattern to identify files that should be renamed
# mapping the mapping of extract ID (names of the vector) to sample code (entries in the vector)
# start position of the first character of the extract ID in the file name
# end position after the last character of the extract ID in the file name
replaceExtractId = function(path, pattern, mapping, start, end) {
	if (!require(stringr)) {
		stop("Could not load required package \"stringr\".")
  	}
  
  lapply(list.files(path=path, pattern=pattern, full.names=TRUE),
  			 function(x, mapping, start, end) {
  			 	 extract = str_sub(x, start=start, end=end)
  			 	 file.rename(from=x, 
  			 	 			 to=str_replace(x,
  			 	 			 pattern=extract, 
  			 	 			 replacement=mapping[extract, 1]))
  			 },
  			 mapping=mapping,
  			 start=start,
  			 end=end)
}

# Read raw TCGA data
# path absolute or relative directory name
# pattern file name pattern of the files to read in 
# what vector giving the data types of the columns to read, will be passed on to "scan"
# sampleType either tumor or normal; if "tumor", only samples with tissue types 1 through 9 will
# be read; if "normal", only samples with tissue types 10 and 11 will be read
# Returns a list with one entry per file. The sample barcode is used as element name
readRawData = function(path, pattern, what, sampleType=c("tumor", "normal")) {
	if (!require(stringr)) {
    stop("Could not load required package \"stringr\".")
  }
  
	sampleType = match.arg(sampleType)

	# Assume tumor by default
	typeCodes = 1:9
	if (sampleType == "normal") {
		typeCodes = 10:11
	}
	
	read = list()
	for (f in list.files(path=path, pattern=pattern, full.names=TRUE)) {
		# Get the start of the sample barcode
		loc = str_locate(f, "TCGA")[1, 1]
		# the full barcode
  	barcode = str_sub(f, start=loc, end=loc+11)
  	# and the sample type
  	sampleType = str_sub(f, start=loc+13, end=loc+14)
  	if (as.integer(sampleType) %in% typeCodes) {
  		read[[barcode]] = scan(f, 
  													 what=what, 
  													 skip=1, 
  													 flush=TRUE, 
  													 strip.white=TRUE, 
  													 fill=TRUE, 
  													 sep="\t",
  													 quiet=TRUE)
  	}
  }
  read
}

# Create matrix from RNAseqV2 input
# read the output of a call to "readRawData"
# type either gene or isoform
# Returns the expression matrix with genes/isoforms in rows and samples in columns
exprsRNAseqV2 = function(read, type=c("gene", "isoform")) {
	if (!require(stringr)) {
    stop("Could not load required package \"stringr\".")
  }
  
  type = match.arg(type)
  
	# Get the gene symbol for each gene
	# Column format = symbol|entrezid, e.g. A4GALT|53947
	allGenes = read[[1]]$gene
	if (type == "gene") {
		allGenes = sapply(read[[1]]$gene, function(x) { str_sub(x, start=1, end=str_locate(x, "\\|")[1,1]-1) } )
	}
	
	indexes = which(allGenes != "?")
	
	# Create the matrix containing the expression values
	exprs = matrix(0.0, nrow=length(allGenes[indexes]), ncol=length(read))
	for (i in 1:length(read)) {
	  exprs[, i] = read[[i]]$normalized_count[indexes]
	}
	
	rownames(exprs) = allGenes[indexes]
	colnames(exprs) = names(read)

	exprs
}	

# Create matrix from RNAseq input
# read the output of a call to "readRawData"
# type either gene or isoform
# measure which gene expression measure to use (typically raw_counts, median_length_normalized, rpkm)
# Returns the expression matrix with genes/isoforms in rows and samples in columns
exprsRNAseq = function(read, measure) {
	if (!require(stringr)) {
		stop("Could not load required package \"stringr\".")
	}
	
	# Get the gene symbol for each gene
	# Column format = symbol|entrezid, e.g. A4GALT|53947
	allGenes = sapply(read[[1]]$gene, function(x) { str_sub(x, start=1, end=str_locate(x, "\\|")[1,1]-1) } )
		
	indexes = which(allGenes != "?")
	
	# Create the matrix containing the expression values
	exprs = matrix(0.0, nrow=length(allGenes[indexes]), ncol=length(read))
	for (i in 1:length(read)) {
		exprs[, i] = read[[i]][[measure]][indexes]
	}
	
	rownames(exprs) = allGenes[indexes]
	colnames(exprs) = names(read)
	
	exprs
}	

# Create matrix from Illumina Infinium array input
# read the output of a call to "readRawData"
# Returns a list with two elements. 
# The first element is a matrix with methylation beta values with 
# probes in rows and samples in columns. The second element is an
# annotation matrix with gene symbol, chromosome and position.
infiniumMeth = function(read) {
	# Create the methylation matrix
	meth = matrix(NA, nrow=length(read[[1]]$genesym), ncol=length(read))
	for (i in 1:length(read))
	  meth[, i] = read[[i]]$beta
	
	# Assign sample names
	colnames(meth) = names(read)
	
	# Annotation of the probes
	meth.probeann = matrix("", nrow=nrow(meth), ncol=3)
	colnames(meth.probeann) = c("genesym", "chrom", "pos")
	meth.probeann[, "genesym"] = read[[1]]$genesym
	meth.probeann[, "chrom"] = read[[1]]$chrom
	meth.probeann[, "pos"] = read[[1]]$pos
	if ("probe" %in% names(read[[1]])) {
		rownames(meth.probeann) = read[[1]]$probe
	} else {
		rownames(meth.probeann) = paste(read[[1]]$chrom, read[[1]]$pos, sep="_")
	}
	
	rownames(meth) = rownames(meth.probeann)
	
	list(methylation.data=meth, methylation.probeann=meth.probeann)
}

# Create matrix from RNAseqV2 input
# read the output of a call to "readRawData"
# Returns the expression matrix with genes/isoforms in rows and samples in columns
miRNAseq = function(read) {
	indexes = seq(1, length(read[[1]]$miRNA_ID), by=2)
	
	# Get the symbol for each miRNA
	allGenes = read[[1]]$miRNA_ID[indexes]
	
	# Create the matrix containing the expression values
	exprs = matrix(NA, nrow=length(allGenes), ncol=length(read))
	for (i in 1:length(read)) {
		exprs[, i] = read[[i]]$rpkm[indexes][match(read[[i]]$miRNA_ID[indexes], allGenes)]
	}
	
	rownames(exprs) = allGenes
	colnames(exprs) = names(read)
	
	exprs
}	


# Create matrix from SNP6 copy number input
# read the output of a call to "readRawData"
# Returns the copy number segment matrix, each row represents one segment.
snp6CNV = function(read) {
	snp6 = data.frame()
	for (i in 1:length(read)) {
	  snp6 = rbind(snp6,
	               data.frame(barcode=names(read)[i],
	                          chrom=read[[i]]$chrom,
	                          start=read[[i]]$start,
	                          stop=read[[i]]$stop,
	                          num.mark=read[[i]]$num.mark,
	                          seg.mean=read[[i]]$seg.mean))
	}
	
	snp6
}

