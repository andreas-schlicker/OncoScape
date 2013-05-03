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

