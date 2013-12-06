##' Tests whether the given gene is expressed in at least the given percentage of samples
##' @param exprs.mat expression matrix with genes in rows and samples in columns
##' @param gene gene to test
##' @param exprs.threshold threshold above which a gene is considered to be expressed; default: 0;
##' only sensible for log-ratio expression data
##' @param cutoff percent samples that need to show expression for the gene to be considered expressed; 
##' default: 0.5, i.e. considered expressed if expression in 50% of samples is > exprs.threshold
##' @return TRUE/FALSE
##' @author Andreas Schlicker
geneExpressed = function(exprs.mat, gene, exprs.threshold=0, cutoff=0.5) {
	(sum(sapply(exprs.mat[gene, , drop=FALSE], function(x) { x > exprs.threshold })) / ncol(exprs.mat)) >= cutoff
}

##' Convenience function to test for a number of genes whether they are expressed or not
##' @param exprs.mat expression matrix with genes in rows and samples in columns
##' @param gene genes to test
##' @param exprs.threshold threshold above which a gene is considered to be expressed; default: 0;
##' only sensible for log-ratio expression data
##' @param cutoff percent samples that need to show expression for the gene to be considered expressed; 
##' @param default 0.5, i.e. considered expressed if expression in 50% of samples is > exprs.threshold
##' @return a boolean vector
##' @author Andreas Schlicker
genesExpressed = function(exprs.mat, genes, exprs.threshold=0, cutoff=0.5) { 
	sapply(genes, function(gene) { geneExpressed(exprs.mat, gene, exprs.threshold, cutoff) })
}

##' Perform paired or unpaired Wilcoxon tests.
##' Tests are only performed on common features.
##' @param mat1 data matrix with features in rows and samples in columns
##' @param mat2 data matrix with features in rows and samples in columns
##' @param paired boolean indicating if a paired test is to performed; default: TRUE
##' To run paired tests, only samples that occur in both matrices are used.
##' For unpaired tests, all samples are used.
##' @return named vector with p-values, NA if no p-value could be calculated
##' @author Andreas Schlicker
doWilcox = function(mat1, mat2, paired=TRUE) {
	# Get the two groups of samples
	if (paired) {
		common.samples = intersect(colnames(mat1), colnames(mat2))
		mat1 = mat1[, common.samples, drop=FALSE] 
		mat2 = mat2[, common.samples, drop=FALSE]
	}
	
	sapply(intersect(rownames(mat1), rownames(mat2)), 
		   function(x) { tryCatch(wilcox.test(mat1[x, ], mat2[x, ], paired=paired, exact=FALSE)$p.value, error = function(e) NA) })
}

##' Perform Bartlett's test on each row of the given matrix.
##' @param inpMat data matrix with features in rows and samples in columns
##' @param groups a factor indicating which samples belong to the different groups
##' @return vector with p-values
##' @author Andreas Schlicker
doBartlett = function(inpMat, groups=NULL) {
  apply(inpMat, 1, function(x) { tryCatch(bartlett.test(x, g=groups)$p.value, error = function(e) NA) })
}

##' Perform Levene's test as implemented in the lawstat package on each row of the given matrix.
##' @param inpMat data matrix with features in rows and samples in columns
##' @param groups a factor indicating which samples belong to the different groups
##' @param location one of "median", "mean", "trim.mean"
##' @return vector with p-values
##' @author Andreas Schlicker
doLevene = function(inpMat, groups, location=c("median", "mean", "trim.mean")) {
	# Load the package for performing Levene's test
	require(lawstat) || stop("Couldn't load the \"lawstat\" package.")
	# Get the right location
	location = match.arg(location)
	
	apply(inpMat, 1, function(x) { tryCatch(levene.test(x, g=groups)$p.value, error = function(e) NA) })
}

##' Calculate prioritization score by summing over the rows.
##' NA values are ignored.
##' @param mat matrix with data types in the columns
##' @return vector with final scores
##' @author Andreas Schlicker
dataTypeScore = function(mat) {
  rowSums(mat, na.rm=TRUE) 
}

greaterThan = function(x, y) { x > y }
smallerThan = function(x, y) { x < y }

##' Creates a helper function to count the elements greater
##' than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the number of elements greater than the cutoff
##' @export
##' @author Andreas Schlicker
gtCutoff = function(cutoff) { function(x) { sum(x > cutoff) } }

##' Creates a helper function to find the percentage of elements
##' that are greater than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the percentage of elements greater than the cutoff (ranging from
##' 0 to 1.
##' @export
##' @author Andreas Schlicker
gtCutoffPercent = function(cutoff) { function(x) { sum(x > cutoff) / length(x) } }

##' Creates a helper function to count the elements less
##' than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the number of elements less than the cutoff
##' @export
##' @author Andreas Schlicker
ltCutoff = function(cutoff) { function(x) { sum(x < cutoff) } }

##' Creates a helper function to find the percentage of elements
##' that are less than the given cutoff.
##' @param cutoff cutoff value to apply
##' @return the percentage of elements less than the cutoff (ranging from
##' 0 to 1.
##' @export
##' @author Andreas Schlicker
ltCutoffPercent = function(cutoff) { function(x) { sum(x < cutoff) / length(x) } }

##' Convert a vector into a matrix with one row.
##' Names of the vector are preserved as colnames.
##' @param vec input vector
##' @return the new matrix
##' @author Andreas Schlicker
matrixFromVector = function(vec) {
	tmp = matrix(vec, nrow=1)
	colnames(tmp) = names(vec)
	tmp
}

##' Count the number of tumors that differ from the mean in normal samples.
##' All features that are not present in the tumors and/or normals matrix are filtered out.
##' @param feature vector with feature IDs that will be contained in the results
##' @param test.features vector with feature IDs that will be tested; default: features
##' @param tumors numeric matrix with features in rows and tumors samples in columns
##' @param normals numeric matrix with features in rows and normal samples in columns
##' @param regulation either "down" or "up" to test for values lower than or greater than the
##' normal mean
##' @param stddev how many standard deviations does a sample have to be away from the mean
##' @param paired boolean; if TRUE, the value of the tumor sample is compared to the value of
##' its paired normal sample; if FALSE, all tumor sample values are compared to the mean over
##' all normals
##' @return named list with two components; "summary" is a matrix with absolute (1st column) 
##' and relative (2nd column) numbers of affected samples; "samples" is a named list with all 
##' samples affected by a change in this feature
##' @author Andreas Schlicker
countAffectedSamples = function(features, test.features=features, tumors, normals, regulation=c("down", "up"), stddev=1, paired=TRUE) {
	if (length(features) > 0) {
		regulation = match.arg(regulation)
		
		# Get the correct comparison function
		# If we want to find genes with greater expression in tumors get the greaterThan function
		# If we want to find genes with lower expression in tumors, get the smallerThan function
		compare = switch(regulation, down=smallerThan, up=greaterThan)
		stddev = switch(regulation, down=stddev*-1, up=stddev)
		
		if (!is.matrix(normals)) {
			normals = matrixFromVector(normals)
			rownames(normals) = features
		}
		
		if (!is.matrix(tumors)) {
			tumors = matrixFromVector(tumors)
			rownames(tumors) = features
		}
		
		# All features that are supposed to be tested with data for tumors and normals 
		common = intersect(rownames(tumors), intersect(test.features, rownames(normals)))
		# All features that we won't test but that should be contained in the results
		missing = setdiff(features, common)
		
		samples = NULL
		affected = c()
		if (length(common) > 0) {
			matched.samples = intersect(colnames(tumors), colnames(normals))
			if (length(matched.samples) == 0) {
				paired = FALSE
				warning("countAffectedSamples: No paired samples found. Performing unpaired analysis!")
			}
			if (paired) {
				# Which tumor samples have a matched normal?
				normFactor = length(matched.samples)
				# All differences 
				deltaMat = tumors[common, matched.samples, drop=FALSE] - normals[common, matched.samples, drop=FALSE]
			} else {
				# Number of samples to normalize with
				normFactor = ncol(tumors)
				# Calculate comparison value across normals
				deltaMat = tumors[common, , drop=FALSE] - apply(normals[common, , drop=FALSE], 1, mean, na.rm=TRUE)
			}
			# Per gene standard deviation of the differences
			deltaSd = apply(normals[common, , drop=FALSE], 1, sd, na.rm=TRUE)
			# Get the names of the samples that are affected
			samples = apply(compare(deltaMat - stddev*deltaSd, 0), 1, function(x) { names(which(x)) })
			# Count affected samples
			affected = unlist(lapply(samples, length))
		}
		if (is.null(samples)) {
			samples = list()
			samples[common] = c()
		}
		
		# Add all missing features and resort
		affected[missing] = 0
		affected = affected[features]
		samples = samples[features]
		names(samples) = features
		
		# This only happens if no feature was to be tested. 
		# The actual value isn't important but has to be != 0
		if (!exists("normFactor")) {
			normFactor = 1
		}
		
		res = list(summary=cbind(absolute=affected, relative=(affected / normFactor)),
				   samples=samples)
	} else {
		res = list()
	}
	
	res
}

##' Compute the union of all elements of the input list.
##' @samples a named list of vectors to get the union of
##' @use a vector of scores; if given, an entry of 1 indicates
##' that this sample list has to be taken into account; default: NULL
##' Both arguments need to use the same names for elements
##' @return the union
##' @author Andreas Schlicker
sampleUnion = function(samples, use=NULL) {
	indexes = names(samples)
	if (!is.null(use)) {
		indexes = names(which(use == 1))
	}
	
	res = c()
	for (i in indexes) {
		res = union(res, samples[[i]])
	}
	
	res
}

##' Imputes missing data using the impute.knn function as implemented 
##' in the "impute" library. If this library is not installed or impute=FALSE, all
##' probes with missing values will be removed.
##' @param data.data matrix with features in rows and samples in columns
##' @param impute boolean indicating whether missing values should be imputed
##' @param no.na threshold giving the number of missing values from which on a
##' feature will be removed. By default, a feature is only removed if its value is
##' missing for all samples.
##' @return cleaned matrix
##' @author Andreas Schlickers
cleanMatrix = function(data.mat, impute=TRUE, no.na=ncol(data.mat)) {
	# How many values are missing for each probe?
	no.nas = apply(data.mat, 1, function(x) { sum(is.na(x)) })
	if (!require(impute) || !impute) {
		print("Could not load library \"impute\". All probes with missing values will be removed")
		# Remove probes with missing values
		exclude = which(no.nas > 0)
		meth.data.imputed = list(data=data.mat[-exclude, ])
	} else {
		# Probes to be excluded
		exclude = which(no.nas > no.na)
		# Impute missing values
		meth.data.imputed = impute.knn(data.mat[-exclude, ])
	}
	
	meth.data.imputed
}

##' Calculates the difference in mean between tumors and normals.
##' Filters out all features not contained in both matrices.
##' @param tumors matrix with tumor data
##' @param normals matrix with normal data
##' @return vector with differences
##' @author Andreas Schlicker
meanDiff = function(tumors, normals) {
	common.features = intersect(rownames(tumors), rownames(normals))
	apply(tumors[common.features, , drop=FALSE], 1, mean, na.rm=TRUE) - apply(normals[common.features, , drop=FALSE], 1, mean, na.rm=TRUE)
}

##' Filters the two input vectors.
##' @param vec1 first vector of IDs
##' @param vec2 second vector of IDs
##' @param restrict vector of IDs that should be retained; default: NULL (retain all)
##' @param paired retain only IDs that are in both vec1 and vec2
##' @return a list with the two filtered vectors
##' @author Andreas Schlicker
doFilter = function(vec1, vec2, restrict=NULL, paired=TRUE) {
	if (!is.null(restrict)) {
		vec1 = intersect(vec1, restrict)
		vec2 = intersect(vec2, restrict)
	}
	
	if (paired) {
		matched.samples = intersect(vec1, vec2)
		vec1 = matched.samples
		vec2 = matched.samples
	}
	
	list(vec1, vec2)
}

##' Return the vector if "gene" is contained in the list and 
##' an empty vector otherwise.
##' @param gene gene to pull out
##' @param inpList named input list
##' @return vector from the list
##' @author Andreas Schlicker
ifPresent = function(gene, inpList) {
	if (gene %in% names(inpList)) { 
		res = inpList[[gene]] 
	} else {
		res = c()
	}
	
	res
}

##' Adds NA values for missing elements and returns the vector sorted by name.
##' @param genes elements that have to be in the resulting vector
##' @param inpVec vector of interest
##' @return sorted vector
##' @author Andreas Schlicker
sortAddMissing = function(genes, inpVec) {
	missing = setdiff(genes, names(inpVec))
	inpVec[missing] = NA
	inpVec[genes]
}

##' Parses a configuration file into command line arguments format.
##' The resulting vector can be passed on to the option parser.
##' @param config absolute filename of the configuration file
##' @return vector formatted as command line arguments
##' @author Andreas Schlicker
parseConfigFile = function(config) {
	require(stringr) || stop("Can't load required package \"stringr\"!")
	
	paste("--", sapply(readLines(config), 
					   function(x) { if (x != "")  { tokens=str_split(x, "=")[[1]]; 
						   			 paste(str_trim(tokens[1]), str_trim(tokens[2]), sep="=") } }), sep="")
}

##' Creates a list of "OptionParserOption" objects for creating an option parser.
##' @param options vector of option names that should be contained in the list. 
##' Option names are identical to long form of options and to names in the final 
##' parsed list.
##' @return list with OptionParserOption objects
##' @author Andreas Schlicker
getOptionList = function(options) {
	require(optparse) || stop("Can't load required package \"optparse\"!")
	
	all = list(config=make_option(c("-c", "--config"), type="character", help="Configuration file"),
		 	   input=make_option(c("-i", "--input"), type="character", help="RData file with precomputed results"),
			   inputstep2=make_option(c("-u", "--inputstep2"), type="character", help="RData file with precomputed results from step2", default=""), 
		 	   output=make_option(c("-o", "--output"), type="character", help="RData output file for saving results"),
			   genes=make_option(c("-g", "--genes"), type="character", help="Gene file", default=""),
			   genecol=make_option(c("-e", "--genecolumn"), type="integer", help="Column with gene IDs in the gene file", default=1),
			   samples=make_option(c("-s", "--samples"), type="character", help="Sample file", default=""),
			   threads=make_option(c("-t", "--threads"), type="integer", help="Number of parallel threads", default=11),
			   cancer=make_option(c("-a", "--cancers"), type="character", help="Comma-separated list of cancer types to test", default="all"),
			   prefix=make_option(c("-p", "--prefix"), type="character", help="Plotting file prefix", default=""),
			   notopgenes=make_option(c("-n", "--notopgenes"), type="integer", help="Number of top genes", default=10),
			   datatype=make_option(c("-d", "--datatype"), type="character", help="Comma-separated list of data types to test (all, exprs, achilles, cgh, meth, sommut)", default="all"),
			   methylation=make_option(c("-m", "--methylation"), type="character", help="Resolution of methylation data (base or region)", default="base"),
			   filtercol=make_option(c("-f", "--filtercol"), type="character", help="Name of the column with the functional mutation score", default="NULL"),
			   filterscore=make_option(c("-r", "--filterscore"), type="numeric", help="Minimum functional score a mutation needs to have", default="1"))
	   
	all[options]
}

##' Parses command line options. Option values containing commas will be split at them.
##' @param options list with OptionParserOption objects
##' @param args character vector with command line options
##' @return named list with all options
##' @author Andreas Schlicker
parseOptions = function(options, args) {
	require(optparse) || stop("Can't load required package \"optparse\"!")
	require(stringr) || stop("Can't load required package \"stringr\"!")
	
	options = parse_args(OptionParser(option_list=options), args=args)
	lapply(options, function(x) { if (str_detect(x, ",")) str_split(configOptions$cancer, ",")[[1]] else x })
}
