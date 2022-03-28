# ==============================================================================
# Script to summarise slim simulation output across the six scenarios
# Author: Luke Lloyd-Jones
# Date started: 25/03/2022
# Date updated: 25/03/2022
# ==============================================================================
library(e1071)
library(reshape)
library(ggwebthemes)
library(gridExtra)

# Source the plotting functions

source("slim_simulation/rscripts/plot_funcs.R")
source("slim_simulation/rscripts/cluster_func.R")

in.fl <- "slim_simulation/hpc_out/"
sets  <- c("0pt01pc", "0pt001pc", "0pt0001pc",
           "0pt01pc_litter", "0pt001pc_litter", "0pt0001pc_litter")

# ---------------------------------------
# Train the SVM of the first 10 data sets
# for each scenario
# ---------------------------------------

nreps <- 5
for (s in sets)
{
  print(paste0("DOING SET ", s))
  for (i in seq(1, nreps))
  {
    fil.s.i <- paste0(in.fl, s, "/pcrelate_kinship_results_", 
                      s, "_iter_", i, ".rds")
    if (file.exists(fil.s.i))
    {
      k.s.i <- readRDS(fil.s.i)
      print(table(k.s.i$kinship))
      if (i == 1)
      {
        k.5.stck <- k.s.i
      } else {
        k.5.stck <- rbind(k.5.stck, k.s.i)
      }
    }
  }
  save(k.5.stck, 
       file = paste0("slim_simulation/svm_stacks/five_stck", s, ".Rda"))
  # Train the models and write out
  dat    <- data.frame(k.5.stck[, c('kin', 'k0', 'k2')], 
                       y = as.factor(k.5.stck$kinship))
  # Linear 
  svmfit.lin.s <- svm(y ~ ., data = dat, 
                      kernel = "linear", 
                      cost = 10, scale = FALSE)
  save(svmfit.lin.s, 
       file = paste0("slim_simulation/svm_stacks/svm_linear", s, ".Rda"))
  
  # # Radial 
  # svmfit.rad.s <- svm(y ~ ., data = dat,
  #                     kernel = "radial", cost = 5, scale = FALSE)
  # save(svmfit.lin.s, 
  #      file = paste0("slim_simulation/svm_stacks/svm_radial", s, ".Rda"))
}


# ----------------------------------------
# Cycle over the rest of these data and
# do the classification with my clustering
# and the SVM
# ----------------------------------------

nreps <- 50
k <- 1
for (s in sets)
{
  print(paste0("DOING SET ", s))
  svm.s <- get(load(paste0("slim_simulation/svm_stacks/svm_linear", s, ".Rda")))
  
  for (i in seq(6, nreps))
  {
    print(i)
    fil.s.i <- paste0(in.fl, s, "/pcrelate_kinship_results_", 
                      s, "_iter_", i, ".rds")
    if (file.exists(fil.s.i))
    {
      k.s.i <- readRDS(fil.s.i)
      print(table(k.s.i$kinship))
      
      # Cluster
      k.s.i$ClstKin     <- clustKs(k.s.i)
      print(table(k.s.i$ClstKin))
      
      # Linear SVM
      dat <- data.frame(k.s.i[, c('kin', 'k0', 'k2')], 
                        y = as.factor(k.s.i$kinship))
      new <- predict(svm.s, dat)
      k.s.i$ClstLSVMKin <- new
      print(table(k.s.i$ClstLSVMKin))
      
      # Compute the performance measures
      
      # Kinship
      
      k.s.i$kinship[k.s.i$kinship == "FTP" |
                    k.s.i$kinship == "GGP" |  
                    k.s.i$kinship == "HSP"] <- "2ndDeg"
      
      # Cluster
      
      df.km  <- k.s.i[, c("ClstKin",     "kinship")]
      df.km$ClstKin[df.km$ClstKin == "3rdDeg"] <- "U"
      res.clst <- computeStats(df.km)
      clst <- rbind(melt(res.clst$cmstats), 
                    data.frame(X1 = "All", X2 = "Micro-F1", 
                               value = res.clst$micro.f1),
                    data.frame(X1 = "All", X2 = "Macro-F1", 
                               value = res.clst$macro.f1))
      clst$Method <- "K-means"
      clst$Rep    <- i
      clst$Sim    <- s
      
      # SVM
      
      df.svm <- k.s.i[, c("ClstLSVMKin", "kinship")]
      df.svm$ClstLSVMKin <- as.character(df.svm$ClstLSVMKin)
      df.svm$ClstLSVMKin[df.svm$ClstLSVMKin == "FTP"] <- "2ndDeg"
      df.svm$ClstLSVMKin[df.svm$ClstLSVMKin == "GGP"] <- "2ndDeg"
      df.svm$ClstLSVMKin[df.svm$ClstLSVMKin == "HSP"] <- "2ndDeg"
      res.svm <- computeStats(df.svm)
      
      svms <- rbind(melt(res.svm$cmstats), 
                    data.frame(X1 = "All", X2 = "Micro-F1", 
                               value = res.svm$micro.f1),
                    data.frame(X1 = "All", X2 = "Macro-F1", 
                               value = res.svm$macro.f1))
      svms$Method <- "SVM"
      svms$Rep    <- i
      svms$Sim    <- s
      
      # ----
      # Bind
      # ----
      
      if (k == 1)
      {
        perf.meas <- rbind(clst, svms)
      } else {
        perf.meas <- rbind(perf.meas, rbind(clst, svms))
      }
        
      k <- k + 1
      
    }
    
  }
  
  save(perf.meas, file = paste0("slim_simulation/out/", 
                                s, "_performance_meas.Rdata"))
}

# -----------------------------------------------
# Try and make a visualisation of the performance
# measure results 
# -----------------------------------------------

head(perf.meas)
perf.meas$Measure <- paste0(gsub("Class: ", "", perf.meas$X1), " - ", 
                                 perf.meas$X2)
# Remove the micro F1
perf.meas <- perf.meas[-grep("Micro", perf.meas$X2), ]
# put the Macro-F1 at the end 
perf.meas$Measure <- as.factor(perf.meas$Measure)
perf.meas$Measure <- factor(perf.meas$Measure, 
                 levels =  c("U - Precision",  "U - Recall", 
                             "U - Specificity",
                             "2ndDeg - Precision",   
                             "2ndDeg - Recall",      
                             "2ndDeg - Specificity",        
                             "FSP - Precision",      "FSP - Recall",         
                             "FSP - Specificity",    "POP - Precision",      
                             "POP - Recall",         "POP - Specificity",    
                             "All - Macro-F1"))

# Change the names 
perf.meas$Sim <- as.factor(perf.meas$Sim)
levels(perf.meas$Sim) <- c("0.0001", "0.0001 - litter", "0.001", 
                           "0.001 - litter",  "0.01", "0.01 - litter")
perf.meas$Sim <- factor(perf.meas$Sim, 
                            levels =  c("0.01", "0.001", "0.0001",         
                                        "0.01 - litter", "0.001 - litter",
                                        "0.0001 - litter"))

# Remove the micro F1
perf.meas$Method[perf.meas$Method == "K-means"] <- "Euclidean-distance"

p1   <- ggplot(perf.meas , aes(x = factor(Measure), y = value, 
                               group = factor(Measure)))
p1.1 <- p1 + geom_boxplot(aes(fill = factor(Measure))) + 
          facet_grid(Sim~Method) +
          theme_web_bw() +
          ylab("") + xlab("") + 
          theme(text  = element_text(size = 25, face = "bold"),
                axis.text.x = element_text(size = 10, face = "bold", angle = 45, 
                                           vjust = 1, hjust=1), 
                axis.text.y = element_text(size = 15, face = "bold"),
                legend.position = "none")

png(filename = paste0("slim_simulation/figures/clust_svm_perf_meas.png"),
    bg = "white", width =  12, height = 16, 
    units = 'in', res = 300)
p1.1 
dev.off()

# -----------------------------------
# Example plots of the classification
# -----------------------------------
i <- 6
k.s.i <- readRDS(fil.s.i)
fil.s.i <- paste0(in.fl, s, "/pcrelate_kinship_results_", 
                  s, "_iter_", i, ".rds")
# Cluster
k.s.i$ClstKin     <- clustKs(k.s.i)
print(table(k.s.i$ClstKin))

# Linear SVM
dat <- data.frame(k.s.i[, c('kin', 'k0', 'k2')], 
                  y = as.factor(k.s.i$kinship))
new <- predict(svm.s, dat)
k.s.i$ClstLSVMKin <- new

# Draw the plots

kin.t  <- drawkink0(k.s.i, clstr = F)
k0k2.t <- drawk2k0(k.s.i,  clstr = F)

k.s.i.2 <- k.s.i
k.s.i.2$kinship <- k.s.i.2$ClstKin

kin.g.t  <- drawkink0(k.s.i.2, clstr = T)
k0k2.g.t <- drawk2k0(k.s.i.2,  clstr = T)

k.s.i.2$kinship <- k.s.i.2$ClstLSVMKin
kin.g.c  <- drawkink0(k.s.i.2, clstr = F)
k0k2.g.c <- drawk2k0(k.s.i.2,  clstr = F)

png(filename = paste0("slim_simulation/figures/example_clustering.png"),
    bg = "white", width =  14, height = 16, 
    units = 'in', res = 300)
grid.arrange(kin.t, k0k2.t,
             kin.g.t, k0k2.g.t,
             kin.g.c, k0k2.g.c, ncol = 2)
dev.off()
