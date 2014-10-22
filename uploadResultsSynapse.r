# Upload prioritization results to Synapse
# 
# Author: Andreas Schlicker (a.schlicker@nki.nl)
###############################################################################

uploadFile = function(file, parent, resources) {
	synFile = File(path=file, parentId=parent)
	
	failed = TRUE
	tries = 0
	while (failed && (tries < 5)) {
		res = tryCatch(synStore(synFile, activity=Activity(used=resources)), error=function(e) { print(e); NA })
		if (class(res) == "File") {
			failed=FALSE
		}
		tries = tries + 1
	}
	
	res
}

library(synapseClient)
library(rGithubClient)

synapseLogin()

# Folder with all public datasets
mainFolder = "syn2518467"
tcgaFolder = "syn2518474"
ccleFolder = "syn2518475"
cgcFolder = "syn2518482"
tsgeneFolder = "syn2518483"

# GitHub repository
setGithubAuth("andreas-schlicker", "")
github = getRepo("andreas-schlicker/OncoScape", ref="branch", refName="version1")

# Files for prioritization
pkgFiles = c("achilles.r", "cnamap.r", "expression.r", 
		"methylationmap.r", "preprocessTCGA.r", 
		"scoring.r", "sommut.r", "utils.r")
pkgResources = list()
for (n in pkgFiles) {
	link = getPermlink(github, n)
	pkgResources = append(pkgResources, list(list(url=link, name=basename(link), wasExecuted=FALSE)))
}

tcgaFiles = c("syn1357159", "syn1358293", "syn1356831", "syn1571514", "syn1358417", "syn1356784", 
		"syn1358109", "syn1571428", "syn1571476", "syn1356625", "syn1356968", "syn1445905", 
		"syn1445932", "syn1445879", "syn1571516", "syn1445942", "syn1445871", "syn1445917", 
		"syn1571430", "syn1571478", "syn1445863", "syn1445890", "syn418033", "syn418082", 
		"syn417855", "syn1571504", "syn1446289", "syn417828", "syn418048", "syn1571420",
		"syn1571468", "syn417812", "syn417925", "syn1446246", "syn1446269", "syn1446205", 
		"syn1571506", "syn1446291", "syn1446190", "syn1446255", "syn1571422", "syn1571470", 
		"syn1446185", "syn1446228", "syn415775", "syn416196", "syn412309", "syn1571509", 
		"syn416206", "syn411993", "syn415945", "syn1571424", "syn1571472", "syn411487", 
		"syn412860", "syn1446002", "syn1446017", "syn1445977", "syn1571511", "syn1446028",
		"syn1445965", "syn1446006", "syn1571426", "syn1571474", "syn1445959", "syn1445987", "syn1710680")
tcgaResources = list()
for (n in tcgaFiles) {
	tcgaResources = append(tcgaResources, list(list(entity=n, wasExecuted=FALSE)))
}

# Script for executing the 1st step
step1Script = getPermlink(github, "TCGA_PAN/SCRIPTS/prioritize_tcga_step1.r")
step1Resource = list(list(url=step1Script, name=basename(step1Script), wasExecuted=TRUE))
# Script for executing the 2nd step
step2Script = getPermlink(github, "TCGA_PAN/SCRIPTS/prioritize_tcga_step2.r")
step2Resource = list(list(url=step2Script, name=basename(step2Script), wasExecuted=TRUE))


##Upload TCGA result files
step1FileTCGA = uploadFile("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step1.rdata",
		tcgaFolder, 
		c(pkgResources, tcgaResources, step1Resource))
step2FileTCGA = uploadFile("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/prioritize_tcga_pancancer_allgenes_step2.rdata",
		tcgaFolder, 
		c(pkgResources, tcgaResources, list(list(entity=step1FileTCGA$properties$id, wasExecuted=FALSE)), step2Resource))

for (f in list.files("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/ALLGENES/20140416/NORMAL/PRIORITIZATION", 
		pattern=".tsv", full.names=TRUE)) {
	uploadFile(f, tcgaFolder, list(list(entity=step2FileTCGA$properties$id, wasExecuted=FALSE)))
}

##Upload CCLE result files
step1FileCCLE = uploadFile("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CCLE/20140416/prioritize_tcga_pancancer_allgenes_step1.rdata",
		ccleFolder, 
		c(pkgResources, tcgaResources, step1Resource))
step2FileCCLE = uploadFile("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CCLE/20140416/prioritize_tcga_pancancer_allgenes_step2.rdata",
		ccleFolder, 
		c(pkgResources, tcgaResources, list(list(entity=step1FileCCLE$properties$id, wasExecuted=FALSE)), step2Resource))

for (f in list.files("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CCLE/20140416/PRIORITIZATION", 
		pattern=".tsv", full.names=TRUE)) {
	uploadFile(f, ccleFolder, list(list(entity=step2FileCCLE$properties$id, wasExecuted=FALSE)))
}


## Upload CGC result file 
step2FileCGC = uploadFile("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CGC/20140416/prioritize_tcga_pancancer_cgc.rdata", 
						  cgcFolder,
						  c(pkgResources, tcgaResources, list(list(entity=step1FileTCGA$properties$id, wasExecuted=FALSE)), step2Resource))

for (f in list.files("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CGC/20140416/PRIORITIZATION", 
		pattern=".tsv", full.names=TRUE)) {
	uploadFile(f, cgcFolder, list(list(entity=step2FileCGC$properties$id, wasExecuted=FALSE)))
}

## Upload TSGene result file
step2FileTSGene = uploadFile("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/TSGENES/20140416/prioritize_tsgene.rdata", 
							 tsgeneFolder,
					 		 c(pkgResources, tcgaResources, list(list(entity=step1FileTCGA$properties$id, wasExecuted=FALSE)), step2Resource)) 

for (f in list.files("/srv/nfs4/medoid-bulk/NKI/a.schlicker/PRIORITIZATION/TCGA_PAN/RESULTS/CGC/20140416/PRIORITIZATION", 
		pattern=".tsv", full.names=TRUE)) {
	uploadFile(f, tsgeneFolder, list(list(entity=step2FileTSGene$properties$id, wasExecuted=FALSE)))
}


synapseLogout()
