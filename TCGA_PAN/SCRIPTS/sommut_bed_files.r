library(stringr)

load("tcga_pancancer4_sommut.rdata")

setwd("/srv/nfs4/void-home/NKI/a.schlicker/funseq-0.1")

all.mutations = list()
tcga.sommut.ids = list()
for (cancType in names(tcga.sommut)) {
	temp = data.frame(CHROM=paste("chr", tcga.sommut[[cancType]][, "Chromosome"], sep=""),
										POS=tcga.sommut[[cancType]][, "Start_Position"]-1,
										END=tcga.sommut[[cancType]][, "End_Position"]-1,
										REF=tcga.sommut[[cancType]][, "Reference_Allele"],
										ALT=tcga.sommut[[cancType]][, "Tumor_Seq_Allele1"],
										ALT2=tcga.sommut[[cancType]][, "Tumor_Seq_Allele2"],
										stringsAsFactors=FALSE)
	
	temp[which(temp[, "REF"] == temp[, "ALT"]), "ALT"] = temp[, "ALT2"]
	temp = cbind(temp, ID=apply(temp[, c("CHROM", "POS", "END", "ALT")], 1, function(x) { paste(str_trim(x), collapse="_") }),
										 QUAL=rep(100, times=nrow(temp)),
										 FILTER=rep("PASS", times=nrow(temp)))
  temp = cbind(temp, INFO=temp[, "ID"])
	tcga.sommut.ids[[cancType]] = temp[, "ID"]
  
	all.mutations[[cancType]] = unique(temp[, c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")])
	rownames(all.mutations[[cancType]]) = all.mutations[[cancType]][, "ID"]
}

wf = file("warnings.txt", "a")
funseq.res = list()
for (cancType in names(all.mutations)) {
	print(cancType)
	last = nrow(all.mutations[[cancType]])
	i = 1
	j = 1
	funseq.res[[cancType]] = vector("list", ceiling(last/4500))
	while (i < last) {
		lastInd = min(i+4500, last)
		tryCatch(
			{ tmpFile = tempfile()
				fileCon = file(tmpFile, "w")
				writeLines(c("##fileformat=VCFv4.0", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"), con=fileCon)
				write.table(all.mutations[[cancType]][i:lastInd, ], file=fileCon, row.names=FALSE, col.names=FALSE, sep="\t", quote=FALSE)
				close(fileCon)
				system(paste("export LD_LIBRARY_PATH=/home/NKI/a.schlicker/VAT/lib:/home/NKI/a.schlicker/LIBS/LIBBIOS/lib:/home/NKI/a.schlicker/LIBS/GD/lib:$LD_LIBRARY_PATH; export PATH=/home/NKI/a.schlicker/VAT/bin:/home/NKI/a.schlicker/LIBS/PERL/bin:$PATH; export PERL5LIB=/home/NKI/a.schlicker/funseq-0.1/FUNSEQ/lib:/home/NKI/a.schlicker/LIBS/PERL/lib/perl5:$PERL5LIB; ./funseq -maf 0 -m 1 -inf vcf -outf vcf -f ", tmpFile, sep=""))
				funseq.res[[cancType]][[j]] = read.table(paste("/home/NKI/a.schlicker/funseq-0.1/out/", basename(tmpFile), ".FunSEQ.vcf", sep=""), sep="\t", quote="", header=FALSE, as.is=TRUE)
				unlink(tmpFile)
				unlink(paste("/home/NKI/a.schlicker/funseq-0.1/out/", basename(tmpFile), ".FunSEQ.vcf", sep=""))
				unlink(paste("/home/NKI/a.schlicker/funseq-0.1/out/", basename(tmpFile), ".err", sep=""))
				unlink(paste("/home/NKI/a.schlicker/funseq-0.1/out/", basename(tmpFile), ".log", sep=""))
			}, error=function(e) { cat(paste(cancType, j, i, lastInd, e, sep=" "), file=wf) })
		i = lastInd
		j = j + 1
	}
}
close(wf)

funseq.res.com = lapply(funseq.res, function(x) { unique(do.call("rbind", x)) })

funseq.res.details = vector("list", length(funseq.res.com))
names(funseq.res.details) = names(funseq.res.com)
for (x in names(funseq.res.com)) {
	funseq.res.details[[x]] = t(apply(funseq.res.com[[x]], 1, function(x) { score = str_split(x[8], "CDSS=")[[1]][2]; st = "CDSS" ; if (is.na(score)) { score = str_split(x[8], "NCDS=")[[1]][2] ; st = "NCDS" } ; c(x[1:7], st, score) }))
	colnames(funseq.res.details[[x]]) = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "SCORETYPE", "SCORE")
	funseq.res.details[[x]] = data.frame(funseq.res.details[[x]], stringsAsFactors=FALSE)
	funseq.res.details[[x]][, "POS"] = as.integer(funseq.res.details[[x]][, "POS"])
	funseq.res.details[[x]][, "QUAL"] = as.integer(funseq.res.details[[x]][, "QUAL"])
	funseq.res.details[[x]][, "SCORE"] = as.integer(funseq.res.details[[x]][, "SCORE"])
}


for (cancType in names(tcga.sommut)) {
	scores = rep(NA, times=nrow(tcga.sommut[[cancType]]))
	for (i in 1:length(scores)) {
		scores[i] = max(funseq.res.details[[cancType]][which(funseq.res.details[[cancType]][, "ID"] == tcga.sommut.ids[[cancType]][i]), "SCORE"])
	}
	scores[is.infinite(scores)] = NA
	tcga.sommut[[cancType]] = cbind(tcga.sommut[[cancType]], FunSeq.Score=scores)
}

save(tcga.sommut, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer4_sommut.rdata")
save(funseq.res.details, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer_funseq_details.rdata")
save(funseq.res, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer_funseq.rdata")
save(tcga.sommut.ids, file="/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RDATA/tcga_pancancer_sommut_ids.rdata")




findMismatch = function(x) {
	dups = funseq.res.details[[x]][duplicated(funseq.res.details[[x]][, "ID"]), "ID"]
	for (d in dups) {
		if (length(unique(funseq.res.details[[x]][which(funseq.res.details[[x]][, "ID"] == d), "SCORE"])) > 1) {
			print(paste(x, d))
			print(funseq.res.details[[x]][which(funseq.res.details[[x]][, "ID"] == d),])
		}
	}
}



> lapply(names(funseq.res.details), findMismatch)

