# Perform Wilcoxon tests on the given matrix.
# If the matchedSamples argument is defined, a paired test is performed. In this
# case, the first length(matchedSamples) number of columns need to contain group
# 1 and then all matched samples in group 2 should be given. Samples in groups 1
# and 2 need to be in matched order.
# If a non-paired test is to be performed, groups should be a named vector of
# 1 and 2 indicating which samples belong to group 1 and 2, respectively.
doWilcox = function(inpMat, matchedSamples=NULL, groups=NULL) {
  # Paired Wilcoxon test?
  paired = !is.null(matchedSamples)
 
  # Get the two groups of samples
  if (paired) {
    group1 = 1:length(matchedSamples) 
    group2 = (length(matchedSamples)+1):ncol(inpMat)
  } else {
    group1 = which(colnames(inpMat) %in% names(groups[groups == 1]))
    group2 = which(colnames(inpMat) %in% names(groups[groups == 2]))
  }
  
  wilcox.p = apply(inpMat, 1, function(x) { tryCatch(wilcox.test(x[group1], x[group2], paired=paired, exact=FALSE)$p.value, error = function(e) NA) })
  
  return(wilcox.p)
}

# Perform Bartlett's test on each row of the given matrix.
# Groups should be a named vector indicating which samples belong to the  
# different groups.
doBartlett = function(inpMat, groups=NULL) {
  bartlett.p = apply(inpMat, 1, function(x) { tryCatch(bartlett.test(x, g=groups)$p.value, error = function(e) NA) })
  
  return(bartlett.p)
}

# Perform Levene's test on each row of the given matrix.
# groups is a factor indicating which samples belong to the different groups.
# location one of "median", "mean", "trim.mean"
doLevene = function(inpMat, groups, location=c("median", "mean", "trim.mean")) {
	# Load the package for performing Levene's test
	require(lawstat)
	# Get the right 
	location = match.arg(location)
	
	levene.p = apply(inpMat, 1, function(x) { tryCatch(levene.test(x, g=groups)$p.value, error = function(e) NA) })
  
	return(levene.p)
}

# Calculate prioritization score by scoring each data type
# x should be a matrix with 10 columns
# sig.cutoff is the significance cut-off used
# direction either "down" or "up" for scoring evidence for the respective direction
# of regulation in tumors 
dataTypeScore = function(mat) {
  apply(mat, 1, sum, na.rm=TRUE) 
}

greaterThan = function(x, y) { x > y }
smallerThan = function(x, y) { x < y }



##DEPRECATED
# Perform Wilcoxon tests on the given matrix.
# If the matchedSamples argument is defined, a paired test is performed. In this
# case, the first length(matchedSamples) number of columns need to contain group
# 1 and then all matched samples in group 2 should be given. Samples in groups 1
# and 2 need to be in matched order.
# If a non-paired test is to be performed, groups should be a named vector of
# 1 and 2 indicating which samples belong to group 1 and 2, respectively.
# olddoWilcox = function(inpMat, matchedSamples=NULL, groups=NULL) {
#   # Paired Wilcoxon test?
#   paired = !is.null(matchedSamples)
#  
#   # p-values from tests
#   wilcox.p = rep(-1, times=nrow(inpMat))
#   if (!is.null(rownames(inpMat)))
#     names(wilcox.p) = rownames(inpMat)
#  
#   # Run paired Wilcoxon tests
#   if (paired) {
#     lms = length(matchedSamples)
#     ncolMat = ncol(inpMat)
#     for (i in 1:nrow(inpMat)) {
#       wilcox.p[i] = tryCatch(wilcox.test(x=inpMat[i, 1:lms], y=inpMat[i, (lms+1):ncolMat], paired=TRUE)$p.value, error = function(e) NA)
#     }
#   } else {
#     group1 = names(groups[groups == 1])
#     group2 = names(groups[groups == 2])
#     for (i in 1:nrow(inpMat)) {
#       wilcox.p[i] = tryCatch(wilcox.test(x=inpMat[i, group1], y=inpMat[i, group2], paired=FALSE)$p.value, error = function(e) NA)
#     }
#   }
#   return(wilcox.p)
# }

##DEPRECATED
# Calculate prioritization score by scoring each single dataset
# x should be a matrix with 10 columns
# categoryScore = function(x) {
#   # Final score
#   res = 0
#   # At least one Illumina 27 is higher methylated in tumors than normals
#   # (significant in paired test) and negatively correlated with expression
#   if (!is.na(x[1]) && x[1] > 0) {
#     res = res + 1
#   }
#   # At least one Illumina 27 is higher methylated in tumors than normals
#   # (significant in unpaired test) and negatively correlated with expression
#   if (!is.na(x[2]) && x[2] > 0) {
#     res = res + 1
#   }
#   # At least one Illumina 450 is higher methylated in tumors than normals
#   # (significant in paired test) and negatively correlated with expression
#   if (!is.na(x[3]) && x[3] > 0) {
#     res = res + 1
#   }
#   # At least one Illumina 450 is higher methylated in tumors than normals
#   # (significant in unpaired test) and negatively correlated with expression
#   if (!is.na(x[4]) && x[4] > 0) {
#     res = res + 1
#   }
#   # Copy number is significantly lower in tumors than in paired normals
#   if (!is.na(x[5]) && x[5] < 0) {
#     res = res + 1
#   }
#   # Gene is significantly associated with poorer disease free survival and
#   # higher expression -> shorter survival
#   if (!is.na(x[7]) && x[7] < 0.05 && !is.na(x[8]) && x[8] < 0) {
#     res = res + 1
#   }
#   # Gene is significantly associated with poorer disease specific survival and
#   # higher expression -> shorter survival
#   if (!is.na(x[9]) && x[9] < 0.05 && !is.na(x[10]) && x[10] < 0) {
#     res = res + 1
#   }
#   return(res)
# }

# # Calculate prioritization score by scoring each data type
# # x should be a matrix with 10 columns
# # sig.cutoff is the significance cut-off used
# # direction either "down" or "up" for scoring evidence for the respective direction
# # of regulation in tumors 
# dataTypeScore = function(x, sig.cutoff=0.05, direction=c("down", "up")) {
#   direction = match.arg(direction)
#  
#   # Get the correct comparison function
#   # If we want to find genes down-regulated in tumor, get the greaterThan function
#   # If we want to find genes up-regulated in tumor, get the smallerThan function
#   compare = switch(direction, down=smallerThan, up=greaterThan)
#   
#   # Final score
#   res = 0
#   # At least one Illumina 27 probe is higher methylated in tumors than normals, 
#   # significant in paired or unpaired test and negatively correlated with expression   
#   # At least one Illumina 450 probe is higher methylated in tumors than normals, 
#   # significant in paired or unpaired test and negatively correlated with expression
#   if ((!is.na(x[1]) && x[1] > 0) || (!is.na(x[2]) && x[2] > 0) ||
#       (!is.na(x[3]) && x[3] > 0) || (!is.na(x[4]) && x[4] > 0)) {
#     res = res + 1
#   }
#   # Copy number is significantly lower in tumors than in paired normals
#   if (!is.na(x[7]) && x[7] == 1) {
#     res = res + 1
#   }
#   # Gene is significantly associated with shorter disease free or disease specific 
#   # survival and lower expression -> shorter survival
#   # Negative coefficients indicate lower expression is associated with shorter survival
#   # Positive coefficients indicate higher expression is associated with shorter survival
#   if ((!is.na(x[8]) && compare(x[8], 0) && !is.na(x[9]) && x[9] < sig.cutoff) || 
#       (!is.na(x[10]) && compare(x[10], 0) && !is.na(x[11]) && x[11] < sig.cutoff)) {
#     res = res + 1
#   }
#   return(res)
# }
