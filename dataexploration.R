library(tidyverse)
library(DataExplorer)
library(ggplot2)
library(GGally)

train <- read.csv("train.csv")
test <- read.csv("test.csv")

cancer <- bind_rows(train, test)

plot_missing(cancer)
plot_missing(train)
plot_missing(test)
boxplot(cancer$PSSM) 
boxplot(cancer$SVM) #have all
boxplot(cancer$SVM)
boxplot(cancer$ANN) #have all
ggpairs(cancer)
original.cor <- cor(cancer[,c(-1,-8,-9)],use = "complete.obs", method = "spearman")

plot_correlation(cancer, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))


##############################
########DATA CLEANING#########
##############################

## Amino acid to 0 and 1 (A and S)
cancer <- cancer %>%
  mutate(Amino.Acid = ifelse(Amino.Acid == "A", 0, 1))

## Stochastic reg imputation for PSSM
pssm.lm <- lm(PSSM ~ SVM + ANN, data=cancer)
pssm.preds <- (predict(pssm.lm, newdata=(cancer %>% filter(is.na(PSSM))))+
                  rnorm(sum(is.na(cancer$PSSM)), 0, sigma(pssm.lm)))
cancer <- cancer %>%
  mutate(PSSM=replace(PSSM, is.na(PSSM), pssm.preds))

#check correlation change
plot_correlation(cancer, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))

## Stochastic reg imputation for Consensus
consens.lm <- lm(Consensus ~ SVM + ANN + PSSM, data=cancer)
consens.preds <- (predict(consens.lm, newdata=(cancer %>% filter(is.na(Consensus))))+
                 rnorm(sum(is.na(cancer$Consensus)), 0, sigma(consens.lm)))
cancer <- cancer %>%
  mutate(Consensus=replace(Consensus, is.na(Consensus), consens.preds))

#Check correlation change
plot_correlation(cancer, type="continuous", 
                 cor_args=list(use="pairwise.complete.obs"))

rm(list=c("pssm.lm", "consens.lm"))

