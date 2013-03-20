# Imports TCGA level 3 methylation data.
# directory is a string giving the directory with data to import
importMethLevel3 = function(directory) {
  if (!require(stringr)) {
    stop("Could not load required package \"stringr\".")
  }
  
  batch = list()
  for (f in list.files(path=directory, pattern=".txt$", full.names=TRUE)) {
    # Find the start of the sample barcode
    start = str_locate_all(f, "TCGA")
    start = start[[1]][nrow(start[[1]]), 1]
    # Read the file and use the sample barcode as index
    batch[[substr(f, start=start, stop=start+14)]] = scan(f, what=list(barcode="", beta=0.0, genesym="", chrom="", pos=0), skip=1, flush=TRUE, strip.white=TRUE, fill=TRUE, sep="\t", quiet=TRUE)
  }
   
  # Create the data matrix
  batch.data = matrix(0.0, nrow=length(batch[[1]]$genesym), ncol=length(batch))
  for (i in 1:length(batch))
    batch.data[, i] = batch[[i]]$beta
  # Assign sample names
  cns = c()
  for (i in 1:length(batch))
    cns = c(cns, batch[[i]]$barcode[1])
  colnames(batch.data) = cns
  
  # Annotation of the probes
  batch.data.probeann = matrix("", nrow=nrow(batch.data), ncol=3)
  colnames(batch.data.probeann) = c("genesym", "chrom", "pos")
  batch.data.probeann[, "genesym"] = batch[[1]]$genesym
  batch.data.probeann[, "chrom"] = batch[[1]]$chrom
  batch.data.probeann[, "pos"] = batch[[1]]$pos
  rownames(batch.data.probeann) = paste(batch.data.probeann[, "chrom"], batch.data.probeann[, "pos"], sep="_")
  rownames(batch.data) = paste(batch.data.probeann[, "chrom"], batch.data.probeann[, "pos"], sep="_")
  
  return(list(methylation.data=batch.data, methylation.probeann=batch.data.probeann))
}

# Imports TCGA level 3 gene expression data obtained through RNAseq.
# directory is a string giving the directory with data to import.
# datatype should be either "rpkm", "raw_counts", or "median_length_normalize"
importRNAseqExpressionLevel3 = function(directory, datatype="rpkm") {
  if (!require(stringr)) {
    stop("Could not load required package \"stringr\".")
  }
  
  tcga.mrna = list()
  for (f in list.files(path=directory, pattern="gene.txt$", full.names=TRUE)) {
    # Find the start of the sample barcode
    start = str_locate_all(f, "TCGA")
    start = start[[1]][nrow(start[[1]]), 1]
    # Read the file and use the sample barcode as index
    tcga.mrna[[substr(f, start=start, stop=start+14)]] = scan(f, what=list(barcode="", gene="", raw_counts=double(), median_length_normalize=double(), rpkm=double()), skip=1, flush=TRUE, strip.white=TRUE, fill=TRUE, sep="\t", quiet=TRUE)
  }
  
  # Get the gene symbol for each gene
  # Column format = symbol|entrezid, e.g. A4GALT|53947
  allGenes = sapply(tcga.mrna[[1]]$gene, function(x) { substr(x, start=1, stop=str_locate(x, "\\|")[1,1]-1) } )
  
  # Create the matrix containing the expression values
  tcga.exprs = matrix(0.0, nrow=length(allGenes[-which(allGenes == "?")]), ncol=length(tcga.mrna))
  for (i in 1:length(tcga.mrna)) {
    if (length(tcga.mrna[[i]][[datatype]][-which(allGenes == "?")]) == nrow(tcga.exprs)) {
      tcga.exprs[, i] = tcga.mrna[[i]][[datatype]][-which(allGenes == "?")]
    }
  }
  
  rownames(tcga.exprs) = allGenes[-which(allGenes == "?")]
  
  # Add sample names as column names
  # Already cut off all parts that do not identify the participant
  cns = c()
  for (i in 1:length(tcga.mrna))
    cns = c(cns, substr(tcga.mrna[[i]]$barcode[1], start=0, stop=12))
  colnames(tcga.exprs) = cns

  return(tcga.exprs)
}

# Imports TCGA level 3 copy number data.
# directory is a string giving the directory with data to import.
importCNVExpressionLevel3 = function(directory) {
  if (!require(stringr)) {
    stop("Could not load required package \"stringr\".")
  }
  
  # Read all file names
  files = list.files(path=directory, pattern="seg.txt$", full.names=TRUE)
  # Find the start of the sample barcode
  start = str_locate_all(files[1], "TCGA")
  start = start[[1]][nrow(start[[1]]), 1]
  # Find samples with two different normal profiles
  # -10 indicates blood and -11 indicates colon normals; keep the normal
  dups = duplicated(substr(files, start=start, stop=start+11), fromLast=TRUE)

  cgh = list()
  for (f in files[!dups]) {
    # Find the start of the sample barcode
    start = str_locate_all(f, "TCGA")
    start = start[[1]][nrow(start[[1]]), 1]
    # Read the file and use the sample barcode as index
    cgh[[substr(f, start=start, stop=start+14)]] = scan(f, what=list(barcode="", chrom=integer(), start=integer(), stop=integer(), num.mark=integer(), seg.mean=double()), skip=1, flush=TRUE, strip.white=TRUE, fill=TRUE, sep="\t", quiet=TRUE)
  }

  # Create a data frame containing the segmentation results in the DNAcopy CBS format
  tcga.cgh = data.frame(barcode=cgh[[1]]$barcode, chrom=cgh[[1]]$chrom, start=cgh[[1]]$start, stop=cgh[[1]]$stop, num.mark=cgh[[1]]$num.mark, seg.mean=cgh[[1]]$seg.mean)
  for (i in 2:length(cgh)) {
    tcga.cgh = rbind(tcga.cgh, data.frame(barcode=cgh[[i]]$barcode, chrom=cgh[[i]]$chrom, start=cgh[[i]]$start, stop=cgh[[i]]$stop, num.mark=cgh[[i]]$num.mark, seg.mean=cgh[[i]]$seg.mean))
  }

  return(tcga.cgh)
}

