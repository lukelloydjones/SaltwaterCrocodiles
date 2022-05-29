# ==============================================================================
# Script to summarise slim simulation output across the six scenarios
# Author: Luke Lloyd-Jones
# Date started: 25/03/2022
# Date updated: 27/05/2022
#  - End of may wanted to see what the results look like if we use the 
#    method proposed in PCRelate to classify kin
# ==============================================================================
library(e1071)
library(reshape)
library(ggwebthemes)
library(gridExtra)
library(knitr)
library(dplyr)

# Source the plotting functions

source("slim_simulation/rscripts/plot_funcs.R")
source("slim_simulation/rscripts/cluster_func.R")
source("slim_simulation/rscripts/perf_meas_funcs.R")

in.fl <- "slim_simulation/hpc_out/"
sets  <- c("0pt01pc", "0pt001pc", "0pt0001pc",
           "0pt01pc_litter", "0pt001pc_litter", "0pt0001pc_litter")

# ---------------------------------------
# Cycle over and compute the means and 
# variances of each of the kinship statisitcs
# ---------------------------------------

nreps <- 50
k <- 1
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
      pops  <- k.s.i[k.s.i$kinship == "POP",  c("kin", "k0", "k2")]
      fsps  <- k.s.i[k.s.i$kinship == "FSP",  c("kin", "k0", "k2")]
      secnd <- k.s.i[k.s.i$kinship == "HSP" |
                     k.s.i$kinship == "GGP" |
                     k.s.i$kinship == "FTP",  c("kin", "k0", "k2")]
      
      pop.mn.i <- colMeans(pops)
      pop.sd.i <- apply(pops, 2, function(x){sd(x)})
      
      fsp.mn.i <- colMeans(fsps)
      fsp.sd.i <- apply(fsps, 2, function(x){sd(x)})
      
      sec.mn.i <- colMeans(secnd)
      sec.sd.i <- apply(secnd, 2, function(x){sd(x)})
      if (k == 1)
      {
        res.mn.sds <- data.frame(POPKinMn = pop.mn.i[1], POPK0Mn = pop.mn.i[2], POPK2Mn = pop.mn.i[3],
                                 POPKinSd = pop.sd.i[1], POPK0Sd = pop.sd.i[2], POPK2Sd = pop.sd.i[3],
                                 FSPKinMn = fsp.mn.i[1], FSPK0Mn = fsp.mn.i[2], FSPK2Mn = fsp.mn.i[3],
                                 FSPKinSd = fsp.sd.i[1], FSPK0Sd = fsp.sd.i[2], FSPK2Sd = fsp.sd.i[3],
                                 SECKinMn = sec.mn.i[1], SECK0Mn = sec.mn.i[2], SECK2Mn = sec.mn.i[3],
                                 SECKinSd = sec.sd.i[1], SECK0Sd = sec.sd.i[2], SECK2Sd = sec.sd.i[3],
                                 REP = i, SIM = s)
      } else {
        res.mn.sds.i <- data.frame(POPKinMn = pop.mn.i[1], POPK0Mn = pop.mn.i[2], POPK2Mn = pop.mn.i[3],
                                   POPKinSd = pop.sd.i[1], POPK0Sd = pop.sd.i[2], POPK2Sd = pop.sd.i[3],
                                   FSPKinMn = fsp.mn.i[1], FSPK0Mn = fsp.mn.i[2], FSPK2Mn = fsp.mn.i[3],
                                   FSPKinSd = fsp.sd.i[1], FSPK0Sd = fsp.sd.i[2], FSPK2Sd = fsp.sd.i[3],
                                   SECKinMn = sec.mn.i[1], SECK0Mn = sec.mn.i[2], SECK2Mn = sec.mn.i[3],
                                   SECKinSd = sec.sd.i[1], SECK0Sd = sec.sd.i[2], SECK2Sd = sec.sd.i[3],
                                   REP = i, SIM = s)
        res.mn.sds <- rbind(res.mn.sds, res.mn.sds.i)
      }
      k <- k + 1
    }
  }
}

dim(res.mn.sds)

means <- aggregate(res.mn.sds[, 1:18], 
                   list(res.mn.sds$SIM), mean)

means <- means[, -grep("Sd", colnames(means))]     
means[, -1] <- round(means[, -1], 3)
rownames(means) <- means$Group.1
means <- means[c("0pt01pc", "0pt001pc", "0pt0001pc",
                 "0pt01pc_litter", "0pt001pc_litter", "0pt0001pc_litter"), ]
kable(means, format = 'latex')

sdss <- means[, -grep("Mn", colnames(means))]     
sdss[, -1] <- round(sdss[, -1], 3)
rownames(sdss) <- sdss$Group.1
sdss <- sdss[c("0pt01pc", "0pt001pc", "0pt0001pc",
               "0pt01pc_litter", "0pt001pc_litter", "0pt0001pc_litter"), ]
kable(sdss, format = 'latex')


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
      
      df.km  <- k.s.i[, c("ClstKin", "kinship")]
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
      
      # ---------------------------------
      # Manic method from PCRelate method
      #  - Tighter
      # ---------------------------------
      
      pop.ko.mx <- 2^(-11/2)
      fsp.ko.mn <- 2^(-7/2)
      
      fst.ks.mx <- 2^-(7/4)
      fst.ks.mn <- 2^-(9/4)
      
      sec.ks.mx <- 2^-(2 + 1/2)
      #sec.ks.mn <- 2^-(13/4)
      
      k.s.i$ClstManicTgt <- "U"
      pops.log <- (k.s.i$kin >= fst.ks.mn) & 
                  (k.s.i$kin  < fst.ks.mx)  & 
                  (k.s.i$k0   < pop.ko.mx)
      fsps.log <- (k.s.i$kin >= fst.ks.mn) & 
                  (k.s.i$kin  < fst.ks.mx)  & 
                  (k.s.i$k0  >= fsp.ko.mn)
      hsps.log <- (k.s.i$kin >= 0.125) & 
                  (k.s.i$kin <  sec.ks.mx)  
      
      k.s.i$ClstManicTgt[pops.log] <- "POP"
      k.s.i$ClstManicTgt[fsps.log] <- "FSP"
      k.s.i$ClstManicTgt[hsps.log] <- "2ndDeg"
      
      df.man <- k.s.i[, c("estkinship.PC", "kinship")]
      df.man$estkinship.PC <- as.character(df.man$estkinship.PC)
      df.man$estkinship.PC[df.man$estkinship.PC == "HSP"] <- "2ndDeg"
      res.man <- computeStats(df.man)
      
      man.1 <- rbind(melt(res.man$cmstats), 
                    data.frame(X1 = "All", X2 = "Micro-F1", 
                               value = res.man$micro.f1),
                    data.frame(X1 = "All", X2 = "Macro-F1", 
                               value = res.man$macro.f1))
      man.1$Method <- "Manich1"
      man.1$Rep    <- i
      man.1$Sim    <- s
      
      df.man.2 <- k.s.i[, c("ClstManicTgt", "kinship")]
      df.man.2$ClstManicTgt <- as.character(df.man.2$ClstManicTgt)
      df.man.2$ClstManicTgt[df.man.2$ClstManicTgt == "HSP"] <- "2ndDeg"
      res.man2 <- computeStats(df.man.2)
      
      man.2 <- rbind(melt(res.man2$cmstats), 
                    data.frame(X1 = "All", X2 = "Micro-F1", 
                               value = res.man2$micro.f1),
                    data.frame(X1 = "All", X2 = "Macro-F1", 
                               value = res.man2$macro.f1))
      man.2$Method <- "Manich2"
      man.2$Rep    <- i
      man.2$Sim    <- s
      
      # ----
      # Bind
      # ----
      
      if (k == 1)
      {
        perf.meas <- rbind(clst, svms, man.1, man.2)
      } else {
        perf.meas <- rbind(perf.meas, rbind(clst, svms, man.1, man.2))
      }
        
      k <- k + 1
      
    }
    
  }
  
  save(perf.meas, file = paste0("slim_simulation/out/", 
                                s, "_performance_meas_2.Rdata"))
}

save(perf.meas, file = paste0("slim_simulation/out/performance_meas_results_2.Rdata"))
# -----------------------------------------------
# Try and make a visualisation of the performance
# measure results 
# -----------------------------------------------

perf.meas <- get(load("slim_simulation/out/perf_meas_results.Rdata"))
perf.meas <- get(load("slim_simulation/out/performance_meas_results_2.Rdata"))
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

levels(perf.meas$Sim) <- c("0.01", "0.001", "0.0001",         
                           "0.01 - clutch", "0.001 - clutch",
                           "0.0001 - clutch")
# Remove the micro F1
perf.meas$Method[perf.meas$Method == "K-means"] <- "Euclidean-distance"
perf.meas$Method[perf.meas$Method == "Manich1"] <- "Manichaikul"
perf.meas <- perf.meas[-which(perf.meas$Method == "Manich2"), ]

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

# Save
save(perf.meas, file = "slim_simulation/out/perf_meas_results.Rdata")

# More summaries

perf.meas[perf.meas$Measure == "All - Macro-F1", ] %>% 
  group_by(Method, Sim) %>% 
  summarise(Mns = mean(value))

perf.meas[perf.meas$Measure == "All - Macro-F1", ] %>% 
  group_by(Method) %>% 
  summarise(Mns = mean(value))

#POP - Precision
perf.meas[perf.meas$Measure == "POP - Precision", ] %>% 
  group_by(Method, Sim) %>% 
  summarise(Mns = mean(value))

perf.meas[perf.meas$Measure == "FSP - Precision", ] %>% 
  group_by(Method, Sim) %>% 
  summarise(Mns = mean(value))

perf.meas[perf.meas$Measure == "2ndDeg - Precision", ] %>% 
  group_by(Method, Sim) %>% 
  summarise(Mns = mean(value))

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
