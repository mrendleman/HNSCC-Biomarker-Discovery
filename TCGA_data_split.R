library('RTCGA')
library('RTCGA.clinical')
library('rlist')
clinpath <- "TCGA-HNSC/gdac.broadinstitute.org_HNSC.Merge_Clinical.Level_1.2016012800.0.0/HNSC.clin.merged.txt"

clinical <- readTCGA(clinpath,'clinical')
survival <- survivalTCGA(clinical)

clinical$patient.new_tumor_events.new_tumor_event.days_to_new_tumor_event_after_initial_treatment
summary(as.factor(clinical$patient.new_tumor_events.new_tumor_event_after_initial_treatment))

## code to identify  variable names in clinical
 library("stringr")
str_inds <- lapply(names(clinical),function(x) {str_detect(x,"hpv")})
vars <- names(clinical)[unlist(str_inds)]
summary(as.factor(clinical[[vars[6]]]))
##


# separate by censor status
surv0 <- survival[which(survival$patient.vital_status==0),]
surv1 <- survival[which(survival$patient.vital_status==1),]

# sort data frames by survival time
surv0 <- surv0[order(surv0$times),]
surv1 <- surv1[order(surv1$times),]

# split into test/train
# function to split data
set.seed(9001)
# returns indices of test data, can be inverted to get train indices
one_third_split <- function(num){
  ind <- c()
  while (length(ind) < num) { # add indices three at a time unless we need to add less than three
    test <- sample(1:3,1) # in each set of three, pick one at random for the test set
    diff <- num-length(ind)
    if (diff >= 3) {
      ind <- append(ind,c(test==1,test==2,test==3))
    } else if (diff ==2){
      ind <- append(ind,c(test==1,test==2))
    }  else {
      ind <- append(ind,c(test==1))
    }
  }
  return(ind)
}


zero_split <- one_third_split(nrow(surv0))
surv0_test <- surv0[zero_split,]
surv0_train <- surv0[!zero_split,]

one_split <- one_third_split(nrow(surv1))
surv1_test <- surv1[one_split,]
surv1_train <- surv1[!one_split,]

data_split = list(train=c(surv0_train$bcr_patient_barcode,surv1_train$bcr_patient_barcode),test=c(surv0_test$bcr_patient_barcode,surv1_test$bcr_patient_barcode))
# save data split to file
saveRDS(data_split,"TCGA_data_split.RDS")
# load data split:
# data_split <- loadRDS("TCGA_data_split.RDA")

# Survival curves for train and test sets
library("survminer")
library("survival")
train_data = rbind(surv0_train,surv1_train)
train_data$grp <- "train"
test_data = rbind(surv0_test,surv1_test)
test_data$grp <- "test"
plots <- list()
trainfit <- survfit(Surv(times, patient.vital_status) ~ 1,
               data = train_data)
plots[[1]] <- ggsurvplot(trainfit, data = train_data, risk.table = TRUE)

testfit <- survfit(Surv(times, patient.vital_status) ~ 1,
                    data = test_data)
plots[[2]] <- ggsurvplot(testfit, data = test_data, risk.table = TRUE)

arrange_ggsurvplots(plots, print = TRUE,ncol = 2, nrow = 1, risk.table.height = 0.4)

full_data <- rbind(train_data,test_data)
fullfit <- survfit(Surv(times, patient.vital_status) ~ grp,
                    data = full_data)
(res.cox <- coxph(Surv(times, patient.vital_status) ~ grp,data = full_data))
ggsurv <- ggsurvplot(fullfit, data = full_data, risk.table = TRUE)
ggsurv$plot <- ggsurv$plot+ 
  ggplot2::annotate("text", 
                    x = 100, y = 0.2, # x and y coordinates of the text
                    label = "Hazard ratio for grouping", size = 5)
ggsurv


