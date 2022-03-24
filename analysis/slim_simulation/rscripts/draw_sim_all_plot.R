# ==============================================================================
# Script to draw the multi-panel scatter plots of the pca, kin/k0, k2/k0 
# plots for both the standard and then multiple litter plots.
# Authors: Luke Lloyd-Jones
# Date started: 24/03/2022
# Date updated: 24/03/2022
# ==============================================================================
library(ggplot2)
library(ggwebthemes)
library(RColorBrewer)

# ---------------
# Plot fucntions
# ---------------

in.fl    <- "slim_simulation/out/"
condition <- "0pt01pc_litter"
ggstruct <- get(load(paste0(in.fl, "pca_struc_", condition, ".rds")))

# ---------------------------
# Plotting functions
# ---------------------------

drawPCA <- function(ggstruct)
{
  g <- ggplot(ggstruct, aes(x = Axis1, y = Axis2, colour = factor(strata))) + 
       geom_point() + 
       theme_web_bw() + scale_colour_brewer(palette="Set2") +
       labs(x = paste0("PC1"), 
            y = paste0("PC2"), 
            color ="Subpopulation") +
       theme(legend.position = "none",
             legend.title = element_text(size = 20, face = "bold"),
             legend.text  = element_text(size = 20, face = "bold"),
             axis.text    = element_text(size = 20, face = "bold"),
             axis.title.x = element_text(size = 22, face = "bold"),
             axis.title.y = element_text(size = 22, face = "bold"))
  return(g)
}

drawkink0 <- function(k)
{
  colrs        <- c("black", "green", "brown", "orange", "gold", "blue")
  names(colrs) <- c("U", "FTP", "GGP", "HSP","FSP","POP")
  k <- k[order(k$kin, decreasing = T), ]
  
  gkink0 <- ggplot(k, aes(x = k0, y = kin, colour = factor(kinship))) +
    geom_segment(aes(x = 0, y = 0, xend = 0, yend = 0.25),
                 linetype = "dashed",
                 col = colrs[names(colrs) == "POP"], lwd = 0.5) +
    geom_segment(aes(x = 0.25, y = 0, xend = 0.25, yend = 0.25),
                 linetype = "dashed",
                 col = colrs[names(colrs) == "FSP"], lwd = 0.5) +
    geom_segment(aes(x = 0, y = 0.25, xend = 0.25, yend = 0.25),
                 linetype = "dashed",
                 col = colrs[names(colrs) == "FSP"], lwd = 0.5) +
    geom_segment(aes(x = 0.5, y = 0, xend = 0.5, yend = 0.125),
                 linetype = "dashed",
                 col = colrs[names(colrs) == "HSP"], lwd = 0.5) +
    geom_segment(aes(x = 0, y = 0.125, xend = 0.5, yend = 0.125),
                 linetype = "dashed",
                 col = colrs[names(colrs) == "HSP"], lwd = 0.5) +
    geom_point(data = subset(k, kinship == 'U'),
               aes(x = k0, y = kin, color = kinship)) +
    geom_point(data = subset(k, kinship == 'POP'),
               aes(x = k0, y = kin, color = kinship)) +
    geom_point(data = subset(k, kinship == 'HSP'),
               aes(x = k0, y = kin, color = kinship)) +
    geom_point(data = subset(k, kinship == 'FSP'),
               aes(x = k0, y = kin, color = kinship)) + 
    geom_point(data = subset(k, kinship == 'FTP'),
               aes(x = k0, y = kin, color = kinship)) + 
    geom_point(data = subset(k, kinship == 'GGP'),
               aes(x = k0, y = kin, color = kinship)) +
    theme_web_bw() + 
    scale_colour_manual(values = colrs) +
    geom_point(aes(x = 0, y = 0.25), 
               col = colrs[names(colrs) == "POP"],  size = 3) +
    geom_point(aes(x = 0.25, y = 0.25), 
               col = colrs[names(colrs) == "FSP"], size = 3) +
    geom_point(aes(x = 0.5, y = 0.125), 
               col = colrs[names(colrs) == "HSP"], size = 3) +
    theme(legend.position = c(0.85, 0.75),
          legend.title = element_text(size = 20, face = "bold"),
          legend.text  = element_text(size = 20, face = "bold"),
          axis.text    = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold")) +
    guides(colour=guide_legend(title = "Kinship")) +
    xlab("k(0)") + 
    ylab("Kinship coefficient")
  return(gkink0)
}

drawk2k0 <- function(k)
{
  colrs        <- c("black", "green", "brown", "orange", "gold", "blue")
  names(colrs) <- c("U", "FTP", "GGP", "HSP","FSP","POP")
  gk0k2 <- ggplot(k, aes(x = k0, y = k2, colour = factor(kinship))) +
    geom_point(aes(x = 0, y = 0), 
               col = colrs[names(colrs) == "POP"],  
               size = 3) +
    geom_point(aes(x = 0.25, y = 0.25), 
               col = colrs[names(colrs) == "FSP"], 
               size = 3) +
    geom_point(aes(x = 0.5,  y = 0), 
               col = colrs[names(colrs) == "HSP"], 
               size = 3) +
    geom_segment(aes(x = 0.25, y = 0, xend = 0.25, yend = 0.25),
                 linetype = "dashed", 
                 col = colrs[names(colrs) == "FSP"], lwd = 0.5) +
    geom_segment(aes(x = 0, y = 0.25, xend = 0.25, yend = 0.25),
                 linetype = "dashed", 
                 col = colrs[names(colrs) == "FSP"], lwd = 0.5) +
    geom_point(data = subset(k, kinship == 'U'),
               aes(x = k0, y = k2, color = kinship)) +
    geom_point(data = subset(k, kinship == 'POP'),
               aes(x = k0, y = k2, color = kinship)) +
    geom_point(data = subset(k, kinship == 'HSP'),
               aes(x = k0, y = k2, color = kinship)) +
    geom_point(data = subset(k, kinship == 'FSP'),
               aes(x = k0, y = k2, color = kinship)) + 
    geom_point(data = subset(k, kinship == 'FTP'),
               aes(x = k0, y = k2, color = kinship)) + 
    geom_point(data = subset(k, kinship == 'GGP'),
               aes(x = k0, y = k2, color = kinship)) +
    theme_web_bw() + 
    scale_colour_manual(values = colrs) +
    theme(legend.position = "none",
          legend.title = element_text(size = 20, face = "bold"),
          legend.text  = element_text(size = 20, face = "bold"),
          axis.text   = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold")) +
    guides(colour=guide_legend(title = "Kinship")) +
    xlab("k(0)") + 
    ylab("k(2)")
  return(gk0k2)
}

# ----------------------
# Standard sim
# ----------------------

# PCA component
condition.1 <- "0pt01pc"
condition.2 <- "0pt001pc"
condition.3 <- "0pt0001pc"
ggstruct.1 <- get(load(paste0(in.fl, "pca_struc_", condition.1, ".rds")))
ggstruct.2 <- get(load(paste0(in.fl, "pca_struc_", condition.2, ".rds")))
ggstruct.3 <- get(load(paste0(in.fl, "pca_struc_", condition.3, ".rds")))
pca.1 <- drawPCA(ggstruct.1)
pca.2 <- drawPCA(ggstruct.2)
pca.3 <- drawPCA(ggstruct.3)

# Kinship components
k.1 <- readRDS(paste0(in.fl, "pcrelate_kinship_results_", 
                              condition.1, ".rds"))
k.2 <- readRDS(paste0(in.fl, "pcrelate_kinship_results_", 
                              condition.2, ".rds"))
k.3 <- readRDS(paste0(in.fl, "pcrelate_kinship_results_", 
                              condition.3, ".rds"))
kk0.1 <- drawkink0(k.1)
kk0.2 <- drawkink0(k.2)
kk0.3 <- drawkink0(k.3)

k2k0.1 <- drawk2k0(k.1)
k2k0.2 <- drawk2k0(k.2)
k2k0.3 <- drawk2k0(k.3)

# ----------------------------------
# Draw the multi panel supp
# plot. Require you to have
# drawn the three PCRelate plots in
# the slim_simulate.R scripts
# ----------------------------------

png(filename = paste0("slim_simulation/figures/all_slim_standard.png"),
    bg = "white", width =  22, height = 18, 
    units = 'in', res = 300)
grid.arrange(pca.1, kk0.1, k2k0.1,
             pca.2, kk0.2, k2k0.2,
             pca.3, kk0.3, k2k0.3, ncol = 3)
dev.off()

# ----------------------
# Litter sim
# ----------------------

# PCA component
condition.1 <- "0pt01pc_litter"
condition.2 <- "0pt001pc_litter"
condition.3 <- "0pt0001pc_litter"
ggstruct.1 <- get(load(paste0(in.fl, "pca_struc_", condition.1, ".rds")))
ggstruct.2 <- get(load(paste0(in.fl, "pca_struc_", condition.2, ".rds")))
ggstruct.3 <- get(load(paste0(in.fl, "pca_struc_", condition.3, ".rds")))
pca.1 <- drawPCA(ggstruct.1)
pca.2 <- drawPCA(ggstruct.2)
pca.3 <- drawPCA(ggstruct.3)

# Kinship components
k.1 <- readRDS(paste0(in.fl, "pcrelate_kinship_results_", 
                      condition.1, ".rds"))
k.2 <- readRDS(paste0(in.fl, "pcrelate_kinship_results_", 
                      condition.2, ".rds"))
k.3 <- readRDS(paste0(in.fl, "pcrelate_kinship_results_", 
                      condition.3, ".rds"))
kk0.1 <- drawkink0(k.1)
kk0.2 <- drawkink0(k.2)
kk0.3 <- drawkink0(k.3)

k2k0.1 <- drawk2k0(k.1)
k2k0.2 <- drawk2k0(k.2)
k2k0.3 <- drawk2k0(k.3)

# ----------------------------------
# Draw the multi panel supp
# plot. Require you to have
# drawn the three PCRelate plots in
# the slim_simulate.R scripts
# ----------------------------------

png(filename = paste0("slim_simulation/figures/all_slim_litter_sim.png"),
    bg = "white", width =  22, height = 18, 
    units = 'in', res = 300)
grid.arrange(pca.1, kk0.1, k2k0.1,
             pca.2, kk0.2, k2k0.2,
             pca.3, kk0.3, k2k0.3, ncol = 3)
dev.off()

