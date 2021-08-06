# TODO: need to finalize file structure for git repo
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("RTCGA")
BiocManager::install("RTCGA.clinical")
library("RTCGA")

# list of most recent available HNSC datasets:
check <- checkTCGA("DataSets","HNSC")

# download clinical data
downloadTCGA(cancerTypes="HNSC",destDir="TCGA-HNSC",date = "2016-01-28")

# How to load clinical data:
# clinpath <- "TCGA-HNSC/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0/HNSC.clin.merged.txt"
# clinical <- readTCGA(clinpath,'clinical')


# download mRNAseq data
downloadTCGA(cancerTypes="HNSC",destDir="TCGA-HNSC",date="2016-01-28",
             dataSet = "HNSC.mRNAseq_Preprocess.Level_3.2016012800.0.0.tar.gz")

# How to load: 
#   - unzip the folder
#   - preprocess
#   - store preprocessed data for loading

mRNA_file = "TCGA-HNSC/gdac.broadinstitute.org_HNSC.mRNAseq_Preprocess.Level_3.2016012800.0.0/HNSC.uncv2.mRNAseq_RSEM_all.txt"
full_data_load <- read.csv(mRNA_file, sep = "\t")
full_data <- data.frame(t(full_data_load[-1]))
colnames(full_data) <- full_data_load[, 1]

# select only rows with primary tumor samples
tumor_data <- full_data[grep("^TCGA.........01",rownames(full_data)),]
# change row names to patient IDs
rownames(tumor_data) <- gsub("^(TCGA........).*","\\1",rownames(tumor_data))
# clean up IDs
rownames(tumor_data) <- gsub("\\.","-",rownames(tumor_data))
saveRDS(tumor_data,file="TCGA_HNSC_tumor_mRNAseq.RDS")

