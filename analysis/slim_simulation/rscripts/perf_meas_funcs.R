# =============================================================
# Some exploration to get at computing performance measures for 
# classifiers
# Author: Luke
# Date: End of March
# =============================================================
library(caret)

# ref.labels <- c(rep("A", 45), rep("B" , 10), 
#                 rep("C", 15), rep("D", 25), 
#                 rep("E", 5))
# predictions <- c(rep("A", 35), rep("E", 5), rep("D", 5),
#                  rep("B", 9), rep("D", 1),
#                  rep("C", 7), rep("B", 5), rep("C", 3),
#                  rep("D", 23), rep("C", 2),
#                  rep("E", 1), rep("A", 2), rep("B", 2))
# df <- data.frame("Prediction" = predictions, 
#                  "Reference"  = ref.labels, 
#                  stringsAsFactors = TRUE)

# ---------
# Functions
# ---------

calculate.accuracy <- function(predictions, ref.labels) 
{
  return(length(which(predictions == ref.labels)) / length(ref.labels))
}

calculate.w.accuracy <- function(predictions, ref.labels, weights) 
{
  lvls <- levels(ref.labels)
  if (length(weights) != length(lvls)) {
    stop("Number of weights should agree with the number of classes.")
  }
  if (sum(weights) != 1) {
    stop("Weights do not sum to 1")
  }
  accs <- lapply(lvls, function(x) {
    idx <- which(ref.labels == x)
    return(calculate.accuracy(predictions[idx], ref.labels[idx]))
  })
  acc <- mean(unlist(accs))
  return(acc)
}

# -----------------------------
# F1 F1 F1 F1 F1 F1 micro/macro
# -----------------------------

get.conf.stats <- function(cm) 
{
  out <- vector("list", length(cm))
  for (i in seq_along(cm)) {
    x <- cm[[i]]
    tp <- x$table[x$positive, x$positive] 
    fp <- sum(x$table[x$positive, colnames(x$table) != x$positive])
    fn <- sum(x$table[colnames(x$table) != x$positie, x$positive])
    # TNs are not well-defined for one-vs-all approach
    elem <- c(tp = tp, fp = fp, fn = fn)
    out[[i]] <- elem
  }
  df <- do.call(rbind, out)
  rownames(df) <- unlist(lapply(cm, function(x) x$positive))
  return(as.data.frame(df))
}

get.micro.f1 <- function(cm) 
{
  cm.summary <- get.conf.stats(cm)
  tp <- sum(cm.summary$tp, na.rm = T)
  fn <- sum(cm.summary$fn, na.rm = T)
  fp <- sum(cm.summary$fp, na.rm = T)
  pr <- tp / (tp + fp)
  re <- tp / (tp + fn)
  f1 <- 2 * ((pr * re) / (pr + re))
  return(f1)
}

# F1 - micro 

get.macro.f1 <- function(cm) 
{
  cs <- cm[[1]]$byClass # a single matrix is sufficient
  cs <- cs[!is.na(cs[, "Recall"]) |
           !is.na(cs[, "Precision"]), ]
  rec <- cs[, "Recall"]
  prs <- cs[, "Precision"]
  re <- sum(rec, na.rm = T) / (nrow(cs))
  pr <- sum(prs, na.rm = T) / (nrow(cs))
  f1 <- 2 * ((re * pr) / (re + pr))
  return(f1)
}

# All stats

computeStats <- function(df)
{
  colnames(df) <- c("Prediction", "Reference")
  
  df$Reference  <- as.factor(df$Reference)
  df$Prediction <- as.factor(df$Prediction)
  
  if (!all(levels(df$Reference) %in% levels(df$Prediction)))
  {
    ext <- levels(df$Reference)[!levels(df$Reference) %in% levels(df$Prediction)]
    levels(df$Prediction) <- c(levels(df$Prediction), ext) 
  }
  
  if (!all(levels(df$Prediction) %in% levels(df$Reference)))
  {
    ext <- levels(df$Prediction)[!levels(df$Prediction) %in% levels(df$Reference)]
    levels(df$Reference) <- c(levels(df$Reference), ext)
  }
  
  cm <- vector("list", length(levels(df$Reference)))
  for (i in seq_along(cm)) 
  {
    positive.class <- levels(df$Reference)[i]
    
    # in the i-th iteration, use the i-th class as the positive class
    cm[[i]] <- confusionMatrix(df$Prediction, 
                               df$Reference, 
                               positive = positive.class)
  }
  
  metrics <- c("Precision", "Recall", "Specificity")
  print(cm[[1]]$byClass[, metrics])
  
  micro.f1 <- get.micro.f1(cm)
  print(paste0("Micro F1 is: ", round(micro.f1, 2)))
  
  macro.f1 <- get.macro.f1(cm)
  print(paste0("Macro F1 is: ", round(macro.f1, 2)))
  
  return(list(cmstats = cm[[1]]$byClass[, metrics], 
              micro.f1 = micro.f1, 
              macro.f1 = macro.f1))
}



# ------------------------------------
# TEST
# Precision, recall and F1-micro/macro
# ------------------------------------

# Accuracy - multi-class accuracy is defined as the average number of 
# correct predictions

# acc <- calculate.accuracy(df$Prediction, df$Reference)
# print(paste0("Accuracy is: ", round(acc, 2)))
# 
# cm <- vector("list", length(levels(df$Reference)))
# for (i in seq_along(cm)) {
#   positive.class <- levels(df$Reference)[i]
#   # in the i-th iteration, use the i-th class as the positive class
#   cm[[i]] <- confusionMatrix(df$Prediction, df$Reference, 
#                              positive = positive.class)
# }
# 
# 
# metrics <- c("Precision", "Recall")
# print(cm[[1]]$byClass[, metrics])
# 
# micro.f1 <- get.micro.f1(cm)
# print(paste0("Micro F1 is: ", round(micro.f1, 2)))
# 
# macro.f1 <- get.macro.f1(cm)
# print(paste0("Macro F1 is: ", round(macro.f1, 2)))
# 
# 
# # ---------------------------------------
# # Let's try for some of the kinship stuff
# # ---------------------------------------
# 
# # -----------------
# # Simple clustering
# # -----------------
# 
# table(k.s.i$kinship)
# k.s.i$kinship[df.km$kinship == "FTP" |
#               df.km$kinship == "GGP" |  
#               df.km$kinship == "HSP"] <- "2ndDeg"
# df.km  <- k.s.i[, c("ClstKin",     "kinship")]
# df.km$ClstKin[df.km$ClstKin == "3rdDeg"] <- "U"
# colnames(df.km) <- c("Prediction", "Reference")
# 
# acc <- calculate.accuracy(df.km$Prediction, df.km$Reference)
# print(paste0("Accuracy is: ", round(acc, 2)))
# 
# df.km$Reference <- as.factor(df.km$Reference)
# df.km$Prediction <- as.factor(df.km$Prediction)
# 
# cm <- vector("list", length(levels(df.km$Reference)))
# for (i in seq_along(cm)) 
# {
#   positive.class <- levels(df.km$Reference)[i]
#   
#   # in the i-th iteration, use the i-th class as the positive class
#   cm[[i]] <- confusionMatrix(df.km$Prediction, 
#                              df.km$Reference, 
#                              positive = positive.class)
# }
# 
# metrics <- c("Precision", "Recall")
# print(cm[[1]]$byClass[, metrics])
# 
# micro.f1 <- get.micro.f1(cm)
# print(paste0("Micro F1 is: ", round(micro.f1, 2)))
# 
# macro.f1 <- get.macro.f1(cm)
# print(paste0("Macro F1 is: ", round(macro.f1, 2)))
# 
# # -----------------
# # SVM clustering
# # -----------------
# 
# df.svm <- k.s.i[, c("ClstLSVMKin", "kinship")]
# df.svm$ClstLSVMKin <- as.character(df.svm$ClstLSVMKin)
# df.svm$ClstLSVMKin[df.svm$ClstLSVMKin == "FTP"] <- "2ndDeg"
# df.svm$ClstLSVMKin[df.svm$ClstLSVMKin == "GGP"] <- "2ndDeg"
# df.svm$ClstLSVMKin[df.svm$ClstLSVMKin == "HSP"] <- "2ndDeg"
# 
# computeStats <- function(df)
# {
#   colnames(df) <- c("Prediction", "Reference")
# 
#   df$Reference  <- as.factor(df$Reference)
#   df$Prediction <- as.factor(df$Prediction)
# 
#   cm <- vector("list", length(levels(df$Reference)))
#   for (i in seq_along(cm)) 
#   {
#   positive.class <- levels(df.svm$Reference)[i]
#   
#   # in the i-th iteration, use the i-th class as the positive class
#   cm[[i]] <- confusionMatrix(df.svm$Prediction, 
#                              df.svm$Reference, 
#                              positive = positive.class)
# }
# 
#   metrics <- c("Precision", "Recall", "Specificity")
#   print(cm[[1]]$byClass[, metrics])
# 
#   micro.f1 <- get.micro.f1(cm)
#   print(paste0("Micro F1 is: ", round(micro.f1, 2)))
# 
#   macro.f1 <- get.macro.f1(cm)
#   print(paste0("Macro F1 is: ", round(macro.f1, 2)))
#   
#   return(list(cmstats = cm[[1]]$byClass[, metrics], 
#               micro.f1 = micro.f1, 
#               macro.f1 = macro.f1))
# }
# 
# computeStats(df.svm)
# --------------
# Extras
# --------------

# # Weighted accuracy - multi-class accuracy is defined as the average number of 
# # correct predictions. Different weights per class done on number of obs
# # per class
# 
# weights <- rep(1 / length(levels(df$Reference)), length(levels(df$Reference)))
# w.acc <- calculate.w.accuracy(df$Prediction, df$Reference, weights)
# print(paste0("Weighted accuracy is: ", round(w.acc, 2)))





