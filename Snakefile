

rule all:
	input:
		"pipeline_final_output_TBD"

rule download_data:
	input:
		"NONE, to handle"
	output:
		"list of data files"
	shell:
		"Rscript insert_download_script_here.R"

rule data_split:
	input:
		"clinical data file"
	output:
		"TCGA_data_split.RDS",
		"test_train_surv_plot.pdf"
	shell:
		"Rscript TCGA_data_split.R"
	

rule coexpression_graph:
	input:
		"expression data",
		"TCGA_data_split.RDA"
	output:
		"graph as a file"
	shell:
		"command goes here"

rule variant_preproc:
	input:
		"mutation data" #note: does NOT care about test/train split
	output:
		"one-hot-encoded variants"
	shell:
		"command goes here"

rule SPCA_construction:
	input:
		"TCGA-HNSC-data/gdac.broadinstitute.org_HNSC.mRNAseq_Preprocess.Level_3.2016012800.0.0/HNSC.uncv2.mRNAseq_RSEM_all.txt"
		"TCGA_data_split.RDS"
	output:
		"TCGA_spca_train.Rda",
		"TCGA_spca_train_sparse.Rda"
	shell:
		"Rscript TCGA_SPCA.R"

rule SPCA_analysis: # step one of SPCA gene selection; identifies important and unimportant SPCs
	input:
		"TCGA_spca_train.Rda",
        "TCGA_spca_train_sparse.Rda",
		"TCGA-HNSC-data/gdac.broadinstitute.org_HNSC.mRNAseq_Preprocess.Level_3.2016012800.0.0/HNSC.uncv2.mRNAseq_RSEM_all.txt",
		"TCGA_data_split.RDS"
	output:
		"important and unimportant SPCs"
	shell:
		"do stuff"

rule SPCA_GOEA:
	input:
		"important and unimportant SPCs"
	output:
		"pathways associated with important SPCs but not unimportant ones"
	shell:
		"run a script probably"

rule DAE_construction: # grid search of hyperparams finds good DAE weights
	input:
		"mrna data",
		"TCGA_data_split.RDA"
	output:
		"DAE_weights file"
	shell:
		"custom script that uses chichi's code"

