library(tidyverse)
library(DataExplorer)
library(ggplot2)
library(GGally)
library(caret)


cancer <- read.csv("cleanedcancer.csv")


#add factors to columns that need it
cancer$Amino.Acid <- as.factor(cancer$Amino.Acid)
cancer$Response <- as.factor(cancer$Response)

#split data into training and test
cancer.train <- cancer %>% filter(!is.na(Response))
cancer.test <- cancer %>% filter(is.na(Response))

#fit svm model on training data

#writing custom function which svm model will maximize
f1 <- function(data, lev = NULL, model = NULL) {
  f1_val <- F1_Score(y_pred = data$pred, y_true = data$obs, positive = lev[1])
  c(F1 = f1_val)
}

svm.model <- train(form=Response~.,
                   data=cancer.train[,! colnames(cancer.train) %in% c("SiteNum") ],
                   method="svmLinear",
                   metric="F1",
                   preProcess = c('center', 'scale'),
                   trControl=trainControl(method="repeatedcv",
                                          number=10, #Number of pieces of your data
                                          repeats=3,
                                          summaryFunction= f1), #repeats=1 = "cv"
                   tuneLength = 10
)

svm.model$bestTune
svm.model$results

#predict on the test data
svm.preds <- data.frame(Id=cancer.test$SiteNum, Predicted=as.logical(
  as.numeric(predict(svm.model, newdata=cancer.test)) - 1))

#output test data for Kaggle.com upload
write_csv(x=svm.preds, path="./svmpredsTeamCK19.csv")