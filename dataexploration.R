library(tidyverse)
library(DataExplorer)
library(ggplot2)
library(GGally)
library(caret)

train <- read.csv("train.csv")
test <- read.csv("test.csv")

#combine training and test set
cancer <- bind_rows(train, test)

#check overall missing values
plot_missing(cancer)

#checking distribution of main measurement statistics
boxplot(cancer$PSSM) #missing data; planning to linearly impute on regression
boxplot(cancer$SVM) #have all
boxplot(cancer$ANN) #have all

#checking correlation
plot_correlation(cancer, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))


###############################
########DATA CLEANING##########
###############################

## Amino acid to 0 and 1 (T and S)
cancer <- cancer %>%
  mutate(Amino.Acid = ifelse(Amino.Acid == "T", 0, 1))

## Stochastic reg imputation for PSSM
pssm.lm <- lm(PSSM ~ SVM + ANN, data=cancer)
pssm.preds <- (predict(pssm.lm, newdata=(cancer %>% filter(is.na(PSSM))))+
                  rnorm(sum(is.na(cancer$PSSM)), 0, sigma(pssm.lm)))
cancer <- cancer %>%
  mutate(PSSM=replace(PSSM, is.na(PSSM), pssm.preds))


#check correlation change
plot_correlation(cancer, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))


#Stochastic reg imputation for Consensus
consens.lm <- lm(Consensus ~ SVM + ANN + PSSM, data=cancer)
consens.preds <- (predict(consens.lm, newdata=(cancer %>% filter(is.na(Consensus))))+
                 rnorm(sum(is.na(cancer$Consensus)), 0, sigma(consens.lm)))
cancer <- cancer %>%
  mutate(Consensus=replace(Consensus, is.na(Consensus), consens.preds))


#Check correlation change
plot_correlation(cancer, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))


#double check missing values have been filled
plot_missing(cancer)


rm(list=c("pssm.lm", "consens.lm"))

#save cleaned data to be put into model
write_csv(x=svm.preds, path="./cleanedcancer.csv")