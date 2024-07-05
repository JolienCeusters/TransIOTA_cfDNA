#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
##  TransIOTA ctDNA - functions  ##
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#

# Packages
library(boot)
library(pROC)
library(SDMTools)
library(boot)
library(logistf)
library(rms)
library(doBy)

#### 1. Sensitivity and specificity ####

Sensitivity <- function(pred, outcome, threshold, data){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  threshold = threshold
  
  Df = data.frame(p = pred, y = outcome, stringsAsFactors = F)
  
  Confusion <- matrix(nrow = 1, ncol = 5)
  Confusion <- data.frame(Confusion)
  colnames(Confusion) <- c('CutOff', 'TN', 'TP', 'FP', 'FN')
  Confusion$CutOff <- threshold
  
  CM <- confusion.matrix(obs = Df$y, pred = Df$p, threshold = threshold)
  Confusion$TN <- CM[1,1]
  Confusion$TP <- CM[2,2]
  Confusion$FP <- CM[2,1]
  Confusion$FN <- CM[1,2]
  
  if(any(Confusion == 0)){
    Confusion$TN     <- CM[1,1] + 0.1
    Confusion$TP     <- CM[2,2] + 0.1
    Confusion$FN     <- CM[1,2] + 0.1
    Confusion$FP     <- CM[2,1] + 0.1
  }
  
  Confusion$Sensitivity  <- Confusion$TP / (Confusion$TP + Confusion$FN)

  # Logit transformation
  Confusion$logit.sens     <- logit(Confusion$Sensitivity)
  Confusion$logit.se.sens  <- sqrt(1/Confusion$TP + 1/Confusion$FN)
  Confusion$logit.LL.sens  <- Confusion$logit.sens - 1.96*Confusion$logit.se.sens
  Confusion$logit.UL.sens  <- Confusion$logit.sens + 1.96*Confusion$logit.se.sens
  
  Confusion$LL.sens <- inv.logit(Confusion$logit.LL.sens)
  Confusion$UL.sens <- inv.logit(Confusion$logit.UL.sens)
  
  return(c(Sensitivity = Confusion$Sensitivity, LL_Sensitivity = Confusion$LL.sens, UL_Sensitivity = Confusion$UL.sens))
}

Specificity <- function(pred, outcome, threshold, data){
  
  arguments <- as.list(match.call())[-1]
  pred = eval(arguments$pred, data)
  outcome <- eval(arguments$outcome, data)
  threshold = threshold
  
  Df = data.frame(p = pred, y = outcome, stringsAsFactors = F)
  
  Confusion <- matrix(nrow = 1, ncol = 5)
  Confusion <- data.frame(Confusion)
  colnames(Confusion) <- c('CutOff', 'TN', 'TP', 'FP', 'FN')
  Confusion$CutOff <- threshold
  
  CM <- confusion.matrix(obs = Df$y, pred = Df$p, threshold = threshold)
  Confusion$TN <- CM[1,1]
  Confusion$TP <- CM[2,2]
  Confusion$FP <- CM[2,1]
  Confusion$FN <- CM[1,2]
  
  if(any(Confusion == 0)){
    Confusion$TN     <- CM[1,1] + 0.1
    Confusion$TP     <- CM[2,2] + 0.1
    Confusion$FN     <- CM[1,2] + 0.1
    Confusion$FP     <- CM[2,1] + 0.1
  }
  
  Confusion$Specificity  <- Confusion$TN / (Confusion$TN + Confusion$FP)

  # Logit transformation
  Confusion$logit.spec     <- logit(Confusion$Specificity)
  Confusion$logit.se.spec  <- sqrt(1/Confusion$TN + 1/Confusion$FP)
  Confusion$logit.LL.spec  <- Confusion$logit.spec - 1.96*Confusion$logit.se.spec
  Confusion$logit.UL.spec  <- Confusion$logit.spec + 1.96*Confusion$logit.se.spec
  
  Confusion$LL.spec <- inv.logit(Confusion$logit.LL.spec)
  Confusion$UL.spec <- inv.logit(Confusion$logit.UL.spec)
  
  return(c(Specificity = Confusion$Specificity, LL_Specificity = Confusion$LL.spec, UL_Specificity = Confusion$UL.spec))
}


#### 2. delta AUC ####

dAUC <- function(data = data, predictor1 = predictor1, predictor2 = predictor2){
  
  # Calculate ROC curve
  roc1 <- roc(data$OutcomeBin, predictor1) # model1$predict
  roc2 <- roc(data$OutcomeBin, predictor2) # model2$predict
  
  # Calculate the difference
  diff <- roc.test(roc2, roc1, paired = TRUE, method = "delong")
  dAUC <- diff$estimate[1] - diff$estimate[2]
  Var  <- ((dAUC)/diff$statistic)^2 
  LL <- diff$conf.int[1]
  UL <- diff$conf.int[2]

  return(c(dAUC = dAUC, Var = Var, LL = LL, UL = UL))
}


#### 3. delta NB ####

diffNB <- function(data, model1, model2, predictor1, predictor2, outcome, threshold, CI = TRUE, N_boot = 2000){
  
  arguments <- as.list(match.call())[-1]
  pred1 = predictor1 
  pred2 = predictor2 
  outcome <- eval(arguments$outcome, data)
  threshold = threshold
  
  Df = data.frame(p1 = pred1, p2 = pred2, y = outcome, stringsAsFactors = F)
  
  Confusion1 <- matrix(nrow = 1, ncol = 5)
  Confusion1 <- data.frame(Confusion1)
  colnames(Confusion1) <- c('CutOff', 'TN', 'TP', 'FP', 'FN')
  Confusion1$CutOff <- threshold
  
  CM1 <- confusion.matrix(obs = Df$y, pred = Df$p1, threshold = threshold)
  Confusion1$TN <- CM1[1,1]
  Confusion1$TP <- CM1[2,2]
  Confusion1$FP <- CM1[2,1]
  Confusion1$FN <- CM1[1,2]
  Confusion1$NB <- Confusion1$TP / nrow(Df) - Confusion1$FP / nrow(Df) * (threshold / (1 - threshold)) # NB = TP / n - FP / n * (threshold / (1 - threshold))
  
  Confusion2 <- matrix(nrow = 1, ncol = 5)
  Confusion2 <- data.frame(Confusion2)
  colnames(Confusion2) <- c('CutOff', 'TN', 'TP', 'FP', 'FN')
  Confusion2$CutOff <- threshold
  
  CM2 <- confusion.matrix(obs = Df$y, pred = Df$p2, threshold = threshold)
  Confusion2$TN <- CM2[1,1]
  Confusion2$TP <- CM2[2,2]
  Confusion2$FP <- CM2[2,1]
  Confusion2$FN <- CM2[1,2]
  Confusion2$NB <- Confusion2$TP / nrow(Df) - Confusion2$FP / nrow(Df) * (threshold / (1 - threshold)) # TP / n - FP / n * (threshold / (1 - threshold))
  
  diffNB <- Confusion2$NB - Confusion1$NB

  diffNB.boot <- matrix(ncol = 1, nrow = N_boot)
  
  ## Calculate CI (with bootstrapping)
  if(CI == TRUE){
    for(i_boot in 1:N_boot){
      set.seed(i_boot + 10000)
      # draw a sample with replacement: ---------------
      m <- sample(x = 1:nrow(data),
                  size = nrow(data),
                  replace = TRUE)
      x_boot <- data[m, ]
      
      ## Fit model
      # Model 1
      M1.boot <- logistf(model1$formula, data = x_boot,
                         control=logistf.control(maxit = 20000))
      # Model 2
      M2.boot <- logistf(model2$formula, data = x_boot,
                         control=logistf.control(maxit = 20000))
      
      predict1  <- M1.boot$predict
      predict2  <- M2.boot$predict
      boot_pred <- cbind(x_boot, predict1, predict2)
      
      ## Calculate NB
      ConfusionB1 <- matrix(nrow = 1, ncol = 5)
      ConfusionB1 <- data.frame(ConfusionB1)
      colnames(ConfusionB1) <- c('CutOff', 'TN', 'TP', 'FP', 'FN')
      ConfusionB1$CutOff <- threshold
      
      CMb1 <- confusion.matrix(obs = boot_pred$OutcomeBin , pred = boot_pred$predict1, threshold = threshold)
      ConfusionB1$TN <- CMb1[1,1]
      ConfusionB1$TP <- CMb1[2,2]
      ConfusionB1$FP <- CMb1[2,1]
      ConfusionB1$FN <- CMb1[1,2]
      ConfusionB1$NB <- ConfusionB1$TP / nrow(boot_pred) - ConfusionB1$FP / nrow(boot_pred) * (threshold / (1 - threshold)) # TP / n - FP / n * (threshold / (1 - threshold))
      
      ConfusionB2 <- matrix(nrow = 1, ncol = 5)
      ConfusionB2 <- data.frame(ConfusionB2)
      colnames(ConfusionB2) <- c('CutOff', 'TN', 'TP', 'FP', 'FN')
      ConfusionB2$CutOff <- threshold
      
      CMb2 <- confusion.matrix(obs = boot_pred$OutcomeBin, pred = boot_pred$predict2, threshold = threshold)
      ConfusionB2$TN <- CMb2[1,1]
      ConfusionB2$TP <- CMb2[2,2]
      ConfusionB2$FP <- CMb2[2,1]
      ConfusionB2$FN <- CMb2[1,2]
      ConfusionB2$NB <- ConfusionB2$TP / nrow(boot_pred) - ConfusionB2$FP / nrow(boot_pred) * (threshold / (1 - threshold)) # TP / n - FP / n * (threshold / (1 - threshold))
      
      ## Calculate difference in NB
      diffNB.boot[i_boot,] <- ConfusionB2$NB - ConfusionB1$NB
    }
    LL_diffNB <- quantile(diffNB.boot[, 1], probs = 0.025) 
    UL_diffNB <- quantile(diffNB.boot[,1], probs = 0.975) 
    
    return(c(dNB = diffNB, LL = LL_diffNB, UL = UL_diffNB))
  }
  
  return(diffNB)
}


#### 4. Optimism ####

Optimism <- function(orig.var, boot.var, Groups){
  
  df = data.frame(orig.var, boot.var, Groups, stringsAsFactors = F)
  df$Optimism <- df$boot.var - df$orig.var
  
  return(summaryBy(Optimism ~ Groups, data = df))
}







