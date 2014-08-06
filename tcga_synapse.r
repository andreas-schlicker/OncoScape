# Couple of functions to get TCGA data out of Synapse.
# The module also defines some parameter lists that are helpful to access
# any TCGA cancer type and many data types. 
#
# Author: schlandi
###############################################################################

##' Builds the query condition from all given parameters.
##' @param params named list with parameters; names have to match annotation
##' elements used in Synapse
##' @return character vector representation of the query condition
##' @author Andreas Schlicker
buildCondition = function(params) {
	paste(unlist(lapply(names(params), function(x) { paste(x, '=="', params[[x]], '"', sep="") })), collapse=" AND ") 
}

##' Gets the Synapse ID associtated with the data set specified by these parameters.
##' @param params named list with parameters; names have to match annotation
##' elements used in Synapse
##' @param parentId the Synapse ID of the parent, disregarded if == NULL; default: NULL
##' @param sampleSetType either "tumor" or "normal"; default: "tumor"
##' @return Synapse ID or NULL if no data set matched the query condition. More than one
##' Synapse ID might match the query.
##' @author Andreas Schlicker
getSynId = function(params, parentId=NULL, sampleSetType=c("tumor", "normal")) {
	if (!is.null(parentId)) {
		params[["parentId"]] = parentId
	}
	
	if (!is.null(params$sampleSetType)) {
		params$sampleSetType = match.arg(sampleSetType)
	}
	
	res = NULL
	tries = 0
	repeat {
		res = tryCatch(synapseQuery(paste('SELECT id FROM entity WHERE ', buildCondition(params), sep="")), 
					   error=function(e) { if (tries >= 5) stop(e) else NA })
		tries = tries + 1
		if (class(res) == "data.frame") {
			break
		}
	}
	
	res[,1]
}

##' Fetch annotation of available data sets of given type. 
##' @param dataType the data type string; get a list of all data types from getDataTypes()
##' @param parentId the Synapse ID of the tumor type entry
##' @param sampleSetType either "tumor" or "normal"
##' @return the annotation of all available data sets as data frame
##' @author Andreas Schlicker
getAnnotation = function(dataType, parentId, sampleSetType=c("tumor", "normal")) {
	params = getDataParams(dataType)
	params[["parentId"]] = parentId
	if (!is.null(params$sampleSetType)) {
		params$sampleSetType = match.arg(sampleSetType)
	}
	
	res = NULL
	tries = 0
	repeat {
		res = tryCatch(synapseQuery(paste('SELECT id, name, sampleSetType, tissueType, platform, freeze, fileType from entity WHERE ', 
						buildCondition(params), 
						sep="")),
					   error=function(e) { if (tries >= 5) stop(e) else NA })
		tries = tries + 1
		if (class(res) == "data.frame") {
			break
		}
	}
	
	res
}

##' Download the data files associated with the given parameters.
##' The function uses synGet, which takes care of caching. Make
##' sure that the cache directory is configured correctly!
##' @param params named list with parameters; names have to match annotation
##' elements used in Synapse
##' @param parentId the Synapse ID of the parent, disregarded if == NULL; default: NULL
##' @param sampleSetType either "tumor" or "normal"; default: "tumor"
##' @return the objects returned by synGet or NULL if no data set matched the query condition
##' @author Andreas Schlicker
getData = function(params, parentId, sampleSetType) {
	synids = getSynId(params, parentId, sampleSetType)
	res = c()
	if (!is.null(synids)) {
		for (synid in synids) {
			tmp = NA
			tries = 0
			repeat {
				tmp = tryCatch(synGet(synid), 
							   error=function(e) { if (tries >= 5) stop(e) else NA })
				tries = tries + 1
				if (!is.na(tmp)) {
					res = c(res, tmp)
					break
				}
			}
		}
	}
	
	res
}

##' Builds the absolute file name from the given object.
##' @param synData object returned by synGet, e.g. through getData
##' @return file names including absolute path to local file
##' @author Andreas Schlicker
getFileName = function(synData) {
	file.path(synData$cacheDir, synData$files)
}

##' Makes a numeric matrix out of the data frame read from 
##' Synapse TCGA data. All rows with rownames == "?" will be removed.
##' @param inpMat input data frame
##' @param rowN Determines how to handle rownames, either "split" or "keep".
##' "split" causes rownames to be split at "|" and only the first element is kept. 
##' This can be used if rownames are formated like "symbol|entrezid". "keep" does
##' nothing to the rownames. If rowN=="split", the "stringr" package needs to be
##' installed. 
##' @param rmDups determines what to do with duplicated rownames. Currently, 
##' only "remove" is implemented, i.e. only the first occurrence will be kept.
##' @return a numeric matrix
##' @author Andreas Schlicker
transformMat = function(inpMat, rowN=c("split", "keep"), rmDups="remove") {
	require(stringr) || stop("Can't load required packages \"stringr\"!")
	
	rowN = match.arg(rowN)
	
	if (rowN == "split" && !require(stringr)) {
		stop("Could not load required package \"stringr\".")
	}
	
	# rmDups == "remove"
	inpMat = inpMat[!duplicated(inpMat[, 1]), ]
	
	# If samples are duplicated, remove the ones that are not primary tumor
	short = str_sub(colnames(inpMat), start=1, end=12)
	if (length(unique(short)) < ncol(inpMat)) {
		names(short) = colnames(inpMat)
		short = short[short %in% short[duplicated(short)]]
		notPrimary = str_sub(names(short), start=14, end=15) != "01"
		inpMat = inpMat[, setdiff(colnames(inpMat), names(short[notPrimary]))]
		
		# In some cases, there are two primary samples in there
		# Remove them
		short = str_sub(colnames(inpMat), start=1, end=12)
		if (length(unique(short)) < ncol(inpMat)) {
			names(short) = colnames(inpMat)
			short = short[short %in% short[duplicated(short)]]
			inpMat = inpMat[, setdiff(colnames(inpMat), names(short))]
		}
	}
	
	# Get the gene symbol for each gene
	# Column format = symbol|entrezid, e.g. A4GALT|53947
	allGenes = inpMat[, 1]
	if (rowN == "split") {
		allGenes = sapply(inpMat[, 1], function(x) { str_sub(x, start=1, end=str_locate(x, "\\|")[1,1]-1) } )
	}
	
	indexes = which(allGenes != "?")
	
	# Create the matrix containing the expression values
	res = data.matrix(inpMat[indexes, 2:ncol(inpMat), drop=FALSE])
	colnames(res) = str_sub(colnames(inpMat)[2:ncol(inpMat)], start=1, end=12)
	rownames(res) = allGenes[indexes]
	
	res
}	

##' Reshapes a mutation input matrix. The resulting matrix has genes in rows and samples in columns.
##' Non-mutated genes have an empty string entry. Mutations are given in the following form:
##' chromosome.name_start_stop_reference_variant_trv.type_amino.acid.change
##' By default all "-" in sample names are replaced with ".". To keep sample names as they are, 
##' set "replace" and "with" to "". 
##' @param inp input matrix as read from the MAF file
##' @param replace string to replace all occurrences of in sample names; default: "-"
##' @param with string all occurrences of "replace" are replaced with; default: "\\.".  
##' @return the mutation matrix
##' @author Andreas Schlicker
transformMutationMat = function(inp, replace="-", with="\\.") {
	genes = unique(inp[, "Hugo_Symbol"])
	samples = str_replace_all(unique(str_sub(inp[, "Tumor_Sample_Barcode"], start=1, end=12)), replace, with)
	muts = matrix("", ncol=length(samples), nrow=length(genes))
	rownames(muts) = genes
	colnames(muts) = samples
	for (i in 1:nrow(inp)) {
		muts[inp[i, "Hugo_Symbol"], str_sub(inp[i, "Tumor_Sample_Barcode"], start=1, end=12)] = 
				paste(inp[i, c("chromosome_name", "start", "stop", "reference", "variant", "trv_type", "amino_acid_change")], collapse="_")
	}
	
	muts
}

##' Fetch the parameter list for a given cancer type.
##' @param cancerType one of the cancer types
##' @return a list with all parameters or NULL if there is no such cancer type.
##' @author Andreas Schlicker
getCancerParams = function(cancerType) {
	tcga.all[[cancerType]]
}

##' Fetch names of all available cancer types.
##' @return acronyms for all known cancer types.
##' @author Andreas Schlicker
getAllCancerTypes = function() {
	names(tcga.all)
} 

##' Fetch names of available core cancer types.
##' @return acronyms for all core cancer types.
##' @author Andreas Schlicker
getCoreCancerTypes = function() {
	names(tcga.core)
}

##' Fetch names of available extra cancer types.
##' @return acronyms for all extra known cancer types.
##' @author Andreas Schlicker
getExtraCancerTypes = function() {
	names(tcga.extra)
}

##' Fetch the parameter list for a given data type.
##' @param dataType one of the defined data types
##' @return a list with all parameters or NULL if there is no such data type.
##' @author Andreas Schlickers
getDataParams = function(dataType) {
	dataTypeParams[[dataType]]
}

##' Fetch all available data types.
##' @return names for all known data types.
##' @author Andreas Schlicker
getDataTypes = function() {
	names(dataTypeParams)
}

##' Get the specified data set(s). If more than one data set is associated
##' with the set of parameters, all of them are returned.
##' @param cancerType either four letter cancer code used by TCGA or a synapse ID
##' @param dataType one of the data types known to getDataParams()
##' @param sampleSet either "tumor" or "nomal"
##' @param rowN Determines how to handle rownames, either "split" or "keep".
##' "split" causes rownames to be split at "|" and only the first element is kept. 
##' This can be used if rownames are formated like "symbol|entrezid". "keep" does
##' nothing to the rownames. 
##' @param addParams named list with additional parameters to restrict the data sets
##' that will be returned; default: list()
##' @return list with matrices for each data file found
##' @author Andreas Schlicker
getMyData = function(cancerType, dataType, sampleSet=c("tumor", "normal"), rowN=c("split", "keep"), addParams=list()) {
	require("stringr") || stop("Could not load required library \"stringr\"!")
	
	rowN = match.arg(rowN)
	
	# Check whether we already have a Synapse ID. Get it if not.
	if (is.na(str_locate(cancerType, "syn")) || str_locate(cancerType, "syn") != 1) {
		cancerType = getSynId(getCancerParams(cancerType))
	}
	# Get the list of parameters for the data type
	dt.params = c(getDataParams(dataType), addParams)
	# Get all files identified by the parameters
	files = getData(dt.params, cancerType, sampleSet)
	# Read all the files
	tmp = lapply(files, 
			function(x) { read.table(getFileName(x),
						sep="\t", quote="", header=TRUE, comment.char="", 
						stringsAsFactors=FALSE, check.names=FALSE) })
	# Transform all matrices and return the resulting list
	lapply(tmp, transformMat, rowN=rowN)
}

##' Get the sommatic mutations for the given cancer type. The returned matrix
##' can be transformed using "transformMutationMat" if needed.
##' @param cancerType four letter cancer code used by TCGA
##' @param source Synapse ID of the combined MAF files; default: "syn1710680"
##' (current as of October 2013
##' @return matrix form of the MAF file
##' @author Andreas Schlicker
getMyMutationData = function(cancerType, source="syn1710680") {
	mut.file = NA
	tries = 0
	repeat {
		mut.file = tryCatch(synGet(source, downloadFile=TRUE), 
				error=function(e) { if (tries >= 5) stop(e) else NA })
		tries = tries + 1
		if (!is.na(mut.file)) {
			break
		}
	}
	
	read.table(file.path(mut.file$cacheDir, 
					     grep(cancerType, mut.file$files, value=TRUE, ignore.case=TRUE)),
			sep="\t", quote="", header=TRUE, comment.char="", 
			stringsAsFactors=FALSE, fill=TRUE)
}

# List with parameters to access any of the TCGA core data sets
# Synapse IDs for any of the data sets can be obtained using the corresponding parameter list. 
tcga.core = list(LUSC=list(acronym="LUSC", parentId="syn300013"), 
				 READ=list(acronym="READ", parentId="syn300013"),
				 GBM=list(acronym="GBM", parentId="syn300013"),
				 BLCA=list(acronym="BLCA", parentId="syn300013"),
				 UCEC=list(acronym="UCEC", parentId="syn300013"),
			  	 COAD=list(acronym="COAD", parentId="syn300013"),
				 OV=list(acronym="OV", parentId="syn300013"),
				 LAML=list(acronym="LAML", parentId="syn300013"),
				 HNSC=list(acronym="HNSC", parentId="syn300013"),
				 LUAD=list(acronym="LUAD", parentId="syn300013"),
				 BRCA=list(acronym="BRCA", parentId="syn300013"),
				 KIRC=list(acronym="KIRC", parentId="syn300013"))

# List with parameters to access any of the TCGA extra data sets
# Synapse IDs for any of the data sets can be obtained using the corresponding parameter list.
tcga.extra = list(CESC=list(acronym="CESC", parentId="syn300013"),
				  PRAD=list(acronym="PRAD", parentId="syn300013"),
		  		  PAAD=list(acronym="PAAD", parentId="syn300013"),
				  THCA=list(acronym="THCA", parentId="syn300013"),
				  STAD=list(acronym="STAD", parentId="syn300013"),
				  KIRP=list(acronym="KIRP", parentId="syn300013"),
				  LGG=list(acronym="LGG", parentId="syn300013"),  
				  LIHC=list(acronym="LIHC", parentId="syn300013"),
				  KICH=list(acronym="KICH", parentId="syn300013"),
				  DLBC=list(acronym="DLBC", parentId="syn300013"),
				  SARC=list(acronym="SARC", parentId="syn300013"),
				  SKCM=list(acronym="SKCM", parentId="syn300013"))

# Convenience list to access any of the TCGA data sets
# Synapse IDs for any of the data sets can be obtained using the corresponding parameter list.
tcga.all = c(tcga.core, tcga.extra)

# Generic parameter lists to obtain specific data sets.
# These should be valid for any tumor type.
dataTypeParams = list(geneexp=list(dataSubType="geneExp", sampleSetType=""),
					  isoformexp=list(dataSubType="isoformExp", sampleSetType=""),
					  cna=list(dataSubType="cna", fileType="genomicMatrix", sampleSetType=""),
					  cnanocnv=list(dataSubType="cna_nocnv", fileType="genomicMatrix", sampleSetType=""),
					  methylation=list(dataSubType="DNAMethylation", sampleSetType=""),
					  mutation=list(dataSubType="mutation"),
					  rppa=list(dataSubType="RPPA", sampleSetType=""),
					  clincal=list(dataSubType="clinical", sampleSetType=""),
					  samplelist=list(dataSubType="whiteList"),
					  mirna=list(dataSubType="miRNA", sampleSetType=""))
			  
##' Download specified data sets. Make sure that the cache directory is properly set before
##' calling this function.
##' @param tissueTypes list with parameter lists defining the tissue types; default: tcga.all
##' @param dataTypes list with parameter lists defining the data types; default: dataTypeParams
##' @param verbose boolean indicating whether to print status updates; download progress will
##' always be shown for each file
downloadData = function(tissueTypes=tcga.core, dataTypes=dataTypeParams, verbose=TRUE) {
	for (tissueType in names(tissueTypes)) {
		if (verbose) {
			print(tissueType)
		}
		synId = getSynId(tissueTypes[[tissueType]])
		for (datatype in names(dataTypes)) {
			if (verbose) {
				print(datatype)
			}
			for (tissue in c("tumor", "normal")) {
				getData(dataTypes[[datatype]], synId, tissue)
			}
		}
	}
}

##' Everything and the kitchen sink for downloading TCGA data from Synapse.
##' Loads the library, logs into Synapse, sets cache directory (if given),
##' downloads all core and extra cancer types and finally logs out again.
##' If .synapseConfig contains login information and cache directory, no
##' input is required. 
##' @param cacheDir directory to which all data will be downloaded. 
##' If set to NULL, the default Synapse directory will be used; default: NULL
##' @author Andreas Schlicker
downloadAll = function(cacheDir=NULL) {
	require(synapseClient) || stop("Can't load the package \"synapseClient\"! Install using source(\"http://depot.sagebase.org/CRAN.R\") ; pkgInstall(\"synapseClient\")")
	
	synapseLogin()
	
	if (!is.null(cacheDir)) {
		synapseCacheDir(cacheDir)
	}
	
	downloadData()
	downloadData(tissueTypes=tcga.extra, dataTypes=lapply(dataTypeParams, function(x) { x$sampleSetType = NULL; x }))
	
	synapseLogout()
}
