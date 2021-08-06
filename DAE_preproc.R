library("dplyr")

# inputs: mrna data RDS file, data split RDS file
# output: csv file for dae
data_split <- readRDS("TCGA_data_split.RDS")
tumor_data <- readRDS("TCGA_HNSC_tumor_mRNAseq.RDS")
# getting just the train data
tumor_train <- tumor_data[rownames(tumor_data) %in% data_split$train,]

# filtering to top 3000 genes by mean absolute deviation to remove highly invariant genes
tumor_MAD <- apply(tumor_train,2,mad)
MAD_frame <- data.frame(tumor_MAD)
MAD_filter <- rownames(slice_max(MAD_frame,order_by=tumor_MAD,n=3000))
tumor_MAD_filter <- tumor_train[,MAD_filter]

saveRDS(MAD_filter,"MAD_gene_filter.RDS")

normalize_col <- function(column) {
  return(column/max(column))
}
tumor_MAD_norm <- apply(tumor_MAD_filter,2,normalize_col)
write.csv(tumor_MAD_norm,"tumor_train_predae.csv")
write.csv(tumor_MAD_filter,"tumor_train_predae_nonorm.csv")

# filter outliers, SCALE and not-scale

outlier_filter <- function(row,threshold) {
  return(sum(row>threshold)==0)
}

data_pts_kept <- function(threshold) {
  return(sum((apply(tumor_MAD_filter,1,function(x){outlier_filter(x,threshold)}))))
}

xvals <- seq(100000,3000000,length.out=100)
data_pt_counts <- lapply(xvals,data_pts_kept)

plot(xvals,data_pt_counts)

outlier_threshold <- 8e5

tumor_MAD_outl_filt <- tumor_MAD_filter[apply(tumor_MAD_filter,1,function(x){outlier_filter(x,outlier_threshold)}),]
write.csv(tumor_MAD_outl_filt,"tumor_train_nonorm_outl.csv")
tumor_outl_norm <- apply(tumor_MAD_outl_filt,2,normalize_col)
write.csv(tumor_outl_norm,"tumor_train_norm_outl.csv")


maxes <- apply(tumor_MAD_filter,2,max)


# After training DAEs over a wide set of parameters, need to select which
# to use based on validation loss
# Prep validation loss data with some vim regex:
#   - :%s/^M/\r/g  (NOTE: ctrl-v ctrl-m gives proper ^M)
#   - :v/100%|##########|/d
#   - remove leading and trailing text
#   - :%s/^/\=printf('%d,', line('.'))
#   - add header

#dae_frame <- read.csv(file="~/dae_params.csv")
dae_frame <- read.csv(file="~/hpchome/dissertation/Denoising-Autoencoder-Analyzer/cmd_line/Model/MAD_Models4/features.csv")
#dae_loss <- read.csv(file="~/dae_final_loss.csv")
dae_loss <- read.csv(file="~/hpchome/dissertation/Denoising-Autoencoder-Analyzer/cmd_line/Model/MAD4_loss.txt")
dae_frame$Valid_Loss <- dae_loss$Training_Loss
dae_frame$Batch_Size <- as.factor(dae_frame$Batch_Size)
dae_frame$Learning_Rate <- as.factor(dae_frame$Learning_Rate)
dae_frame$Corruption_Rate <- as.factor(dae_frame$Corruption_Rate)
dae_frame$Epochs <- as.factor(dae_frame$Epochs)
dae_frame$Encoder_Activation_Function <- as.factor(dae_frame$Encoder_Activation_Function)

dae_cols <- names(dae_frame)[2:6]

dim_sizes <- unlist(lapply(dae_cols,function(x){return(length(levels(dae_frame[[x]])))}))
dae_names <- lapply(dae_cols,function(x){return(levels(dae_frame[[x]]))})

dae_array <- array(dim=dim_sizes,dimnames=dae_names)
for(batch in dae_names[[1]]) {
  for(learn in dae_names[[2]]) {
    for(corrupt in dae_names[[3]]) {
      for(epoch in dae_names[[4]]) {
        for(encode in dae_names[[5]]) {
          ind <- dae_frame$Batch_Size==batch & dae_frame$Learning_Rate==learn & dae_frame$Corruption_Rate==corrupt & dae_frame$Epoch_Size==epoch & dae_frame$Encoder_Activation_Function==encode
          dae_array[batch,learn,corrupt,epoch,encode] <- dae_frame$Valid_Loss[which(ind)]
        }
      }
    }
  }
}

dae_mat <- as.array(dae_frame)
dae_sigm <- dae_array[,,,,"sigmoid"]

sigm_ind <- which(dae_frame$Encoder_Activation_Function=="sigmoid")
dae_sframe <- data.frame(dae_frame$Batch_Size[sigm_ind])
colnames(dae_sframe) <- "Batch_Size"
dae_sframe$Learning_Rate <- dae_frame$Learning_Rate[sigm_ind]
dae_sframe$Corruption_Rate <- dae_frame$Corruption_Rate[sigm_ind]
dae_sframe$Epoch_Size <- dae_frame$Epoch_Size[sigm_ind]
dae_sframe$Valid_Loss <- dae_frame$Valid_Loss[sigm_ind]



library("ggplot2")
library("gridExtra")
plot_batch <- function(x) {
  batch_ind <- which(dae_frame$Batch_Size==x)
  batch_frame <- data.frame(dae_frame$Learning_Rate[batch_ind])
  colnames(batch_frame) <- "Learning_Rate"
  batch_frame$Corruption_Rate <- dae_frame$Corruption_Rate[batch_ind]
  batch_frame$Epoch_Size <- dae_frame$Epochs[batch_ind]
  batch_frame$Valid_Loss <- dae_frame$Valid_Loss[batch_ind]
  
  return(ggplot(batch_frame,mapping=aes(x=Epoch_Size,y=Valid_Loss,color=Corruption_Rate,shape=Learning_Rate)) + 
    geom_point(size=3) + labs(title=paste0("Batch Size: ",x),x="Epoch Size", y="Validation Loss") +
    theme(plot.title=element_text(hjust=0.5,size=14)) + scale_shape_manual(values=c(19,17,15,3,7,8,10,25)) )
}
# dim 1: learning rate: shape
# dim 2: corruption ratio: color
# dim 3: epoch size: x-axis
# values: y-axis

batch_plots <- lapply(levels(dae_frame$Batch_Size),plot_batch)
grid.arrange(grobs=batch_plots,nrow=2,ncol=2)



## Tumor Data Transformation

# filter ALL tumor data down with MAD filter
tumor_MAD_all <- as.matrix(tumor_data[,MAD_filter])

# normalize tumor data based on training data column max values
normalize_col_train <- function(column) {
  return(tumor_MAD_all[,column]/max(tumor_MAD_filter[,column]))
}
tum_MAD_norm <- sapply(1:ncol(tumor_MAD_all),normalize_col_train)

# after selection of DAE models to use:
#   import encoder W and b
#     to export from numpy:
#       b: python numpy_to_csv.py --input ~/.yadlt/data/dae_MAD/dae_model1/dae_model1-enc_b.npy --name dae_mad_1_encb --output dae_mad_1_encb.csv
#       W: python numpy_to_csv.py --input ~/.yadlt/data/dae_MAD/dae_model1/dae_model1-enc_w.npy --name dae_mad_1_encw --output dae_mad_1_encw.csv --mult_lines True
#     to import to R:

weights_path <- "~/hpchome/dissertation/Denoising-Autoencoder-Analyzer/cmd_line/Model/MAD4_weights/"
W_files <- c("dae_MAD4_92_encw.csv","dae_MAD4_292_encw.csv","dae_MAD4_543_encw.csv","dae_MAD4_743_encw.csv")
b_files <- c("dae_MAD4_92_encb.csv","dae_MAD4_292_encb.csv","dae_MAD4_543_encb.csv","dae_MAD4_743_encb.csv")

library("sigmoid")
load_and_transform <- function(i) {
  dae_W <- t(as.matrix(read.csv(file=paste0(weights_path,W_files[i]),header=FALSE)))
  dae_W <- dae_W[-1,]
  dae_W <- apply(dae_W,2,as.numeric)
  dae_b <- t(as.matrix(read.csv(file=paste0(weights_path,b_files[i]),header=FALSE)))
  dae_b <- as.numeric(dae_b[-1])
  
  presigm <- tum_MAD_norm%*%t(dae_W)+dae_b
  dae_transformed = sigmoid(presigm)
  colnames(dae_transformed) <- unlist(lapply(1:100,function(x){paste0("model_",i,"_node_",x)}))
  return(dae_transformed)
}

daes_list <- lapply(1:4,load_and_transform)
full_transformed <- cbind(daes_list[[1]],daes_list[[2]],daes_list[[3]],daes_list[[4]])

# save dae-transformed data

saveRDS(full_transformed,file="dae_nodes_full.RDS")








