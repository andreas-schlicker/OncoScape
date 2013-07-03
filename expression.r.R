# Functions for performing differential expression analysis using edgeR on
# sequencing data. 
#
# Author: Andreas Schlicker
###############################################################################

##' Creates a DGEList object used by other functions.
##' Rownames of the count.table matrix will be used as rownames for the DGEList object
##' @param count.table numeric read count matrix
##' @param feature.ann annotation for features; default: NULL
##' @param remove.zeros boolean, whether features with only zero counts should
##' be removed
##' @return the DGEList object
##' @author Andreas Schlicker
createDGEList = function(count.table, feature.ann=NULL, remove.zeros=FALSE) {
	if (!require(edgeR)) {
		stop("Could not load required package \"edgeR\"!")
	}
	dgel = DGEList(counts=count.table, genes=gene.ann, remove.zeros=remove.zeros)
	rownames(dgel) = rownames(count.table)
	
	dgel
}

##' Filters out features that are not expressed in given number of samples
##' @param dgel the DGEList object
##' @param count.cutoff count per million cutoff value that is applied to each sample
##' @param sample.cutoff at least sample.cutoff many samples need to meet the 
##' count.cutoff for a feature to be kept; default: 0.1 (= 10% of samples)
##' @param relative boolean, whether sample.cutoff is to be interpreted as a relative
##' fraction of samples or an absolute number; default: TRUE
##' @return filtered DGEList object
##' @author Andreas Schlicker
filterDGEList = function(dgel, count.cutoff, sample.cutoff=0.1, relative=TRUE) {
	# Convert absolute numer of samples to relative quantity
	if (!relative) {
		sample.cutoff = sample.cutoff / ncol(dgel)
	}
		
	keep = (rowSums(cpm(dgel) > count.cutoff) / ncol(dgel)) >= sample.cutoff
	dgel = dgel[keep, ]
	# Correct library size after filtering
	dgel$samples$lib.size <- colSums(dgel$counts)
	
	dgel
}

##' Creates a design matrix.
##' @param sample.names character vector with the sample names; will be used as
##' rownames for the design matrix. Has to contain all tumor samples and then all
##' normal samples. Each entry in this vector has to be unique.
##' @param tumors character vector with the names of the tumor samples 
##' @param normals character vector with the names of the normal samples
##' @only.paired boolean, indicates whether only samples with tumor and normal
##' should be included (default: TRUE). If TRUE, normal and tumor of the 
##' same sample need to have the same name
##' @return the design matrix
##' @author Andreas Schlicker
designMatrix = function(sample.names, tumors, normals, only.paired=TRUE) {
	patients = tumors
	if (only.paired) {
		patients = intersect(tumors, normals)
	}
	tissue = factor(c(rep("T", times=length(patients)), 
					  rep("N", times=length(normals))))
	# Add the normal samples
	patients = factor(c(patients, normals))
	design = model.matrix(~patients+tissue)
	rownames(design) = sample.names
	
	design
}

##' Runs the different dispertion estimation steps.
##' @param dgel DGEList object
##' @param design the design matrix
##' @return the updated DGEList object
##' @author Andreas Schlicker
estimateDispertion = function(dgel, design) { 
	# Estimate the common dispersion
	dgel =  estimateGLMCommonDisp(dgel, design, verbose=FALSE)
	# gene-wise dispersion values
	dgel <- estimateGLMTrendedDisp(dgel, design)
	dgel <- estimateGLMTagwiseDisp(dgel, design)
	dgel
}

##' Runs the differential expression analysis.
##' @param dgel DGEList object
##' @param design design matrix
##' @return data frame with the results of the analysis
##' @author Andreas Schlicker
diffExpr = function(dgel, design) {
	fit = glmFit(dgel, design)
	
	topTags(lrt, n=nrow(lrt))[[1]]
}

rawdata <- read.delim("TableS1.txt", check.names=FALSE, stringsAsFactors=FALSE)
# Create the DGEList object
y <- DGEList(counts=rawdata[,4:9], genes=rawdata[,1:3])


# Correct library size after gene selection
y$samples$lib.size <- colSums(y$counts)
# Set the rownames
rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene
y$genes$EntrezGene <- NULL
# Calculate the normalization factor
y <- calcNormFactors(y)
# Check them out
#y$samples

# Build the design matrix
Patient <- factor(c(8,8,33,33,51,51))
Tissue <- factor(c("N","T","N","T","N","T"))
data.frame(Sample=colnames(y),Patient,Tissue)
design <- model.matrix(~Patient+Tissue)
rownames(design) <- colnames(y)

# Estimate the dispersion
y <- estimateGLMCommonDisp(y, design, verbose=TRUE)
# gene-wise dispersion values
y <- estimateGLMTrendedDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)

# Differential experssion
fit <- glmFit(y, design)

# Get the data.frame with the results
head(topTags(lrt, n=nrow(lrt))[[1]])
