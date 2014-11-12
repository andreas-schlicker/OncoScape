OncoScape
=========

OncoScape is a package for gene prioritization in the R statistical programming environment. The analysis is run in a contrast fashion, i.e. always two groups of samples are compared with each other. Examples include:

  - tumors vs. normals
  - cell lines vs. normals
  - treatment responders vs resistant
  - samples with mutations in gene X vs wild type

Currently, analyses of five data types are implemented in OncoScape:
  1. gene expression
  2. DNA copy number
  3. DNA methylation
  4. mutation
  5. shRNA knock-down data
  
Aberrations in each gene are called for each data type seperately and scored as 0 (no aberration found) or 1 (aberration found). These scores are summed across data types to give the final score. OncoScape differentiates between activating (oncogene-like) and inactivating (tumor suppressor-like) aberrations and calculates independent scores for both directions. It is possible to run the analysis on any combination of these data types. 


Analysis workflow
-----------------

The analysis proceeds in four steps:
  0. data import
  1. calculate statistics (using "do*Analysis" methods)
  2. filtering and scoring (using "summarize*" methods)
  3. plotting of results
  
Example scripts for putting everything together are:
  - TCGA_PAN/SCRIPTS/load_tcga_pancancer_data.r
  - TCGA_PAN/SCRIPTS/prioritize_tcga_step1.r
  - TCGA_PAN/SCRIPTS/prioritize_tcga_step2.r
  - TCGA_PAN/SCRIPTS/plot*
  
  
Source files
------------
  - achilles.r: functions to analyse data from Project Achilles (http://www.broadinstitute.org/achilles)
  - cnamap.r: functions to analyse copy number data
  - compare.r: functions for comparing OncoMap results from two different runs
  - expression.r: functions to analyse gene expression data
  - methylationmap.r: functions to analyse DNA methylation data
  - plotting.r: functions for plotting results
  - preprocessTCGA.r: functions for preprocessing data from the TCGA pancancer project
  - scoring.r: functions for calculating final scores
  - sommut.r: functions to analyse somatic mutation data
  - survival.r: functions to analyse survival data (not used at the moment)
  - tcga_synapse.r: functions to download TCGA data from Sage Bionetwork's Synapse (www.synapse.org)
  - uploadResutlsSynapse.r: functions to upload result filees to Synapse
  - utils.r: helper functions used by others in this project
