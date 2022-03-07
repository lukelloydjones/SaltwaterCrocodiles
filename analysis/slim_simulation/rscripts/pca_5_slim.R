# =========================================================
# Script to do a variety of PCA plots for a slim dataset
# Authors: Kira Villiers, adapted from Luke Lloyd-Jones
# Date started: 19/01/2022
# Date updated: 19/01/2022
# =========================================================

library(ade4)
library(adegenet)
library(SNPRelate) #has snpgdsOpen
library(GWASTools) #for gdsSubset
library(data.table)
library(ggplot2)
library(ggwebthemes)
library(RColorBrewer) #display.brewer.all(colorblindFriendly=TRUE)

#"~/Desktop/slim_mg_2p_migr/"
data.path <- "~/Desktop/slim/slim_sim_2/slim_mg_2p_migr_0pt01pc/" 
data.path <- "~/Desktop/slim/slim_sim_2/slim_mg_2p_migr_0pt0001pc/" 

# ------ Paths to be used --------------
showfile.gds(closeall = TRUE)

snpgdsVCF2GDS(vcf.fn = paste0(data.path, "nonWF_st1_gen_6000.vcf"),
              out.fn = paste0(data.path, "nonWF_st1_gen_6000.gds"))

showfile.gds(closeall = TRUE)

f <- snpgdsOpen(paste0(data.path, "nonWF_st1_gen_6000.gds"))
snpset   <- read.gdsn(index.gdsn(f, "snp.id"))
indset   <- read.gdsn(index.gdsn(f, "sample.id"))
geno.mat <- snpgdsGetGeno(f, snpfirstdim = F)
rownames(geno.mat) <- indset
if (any(is.numeric(snpset))) {
  colnames(geno.mat) <- paste0("m",snpset)
} else {
  colnames(geno.mat) <- snpset
}


selectedind <- indset #sample(indset, 200)
SNPrate <- snpgdsSNPRateFreq(f, with.snp.id = TRUE,
                             sample.id = selectedind)
#hist(SNPrate$MinorFreq)
segregating <- SNPrate$snp.id[SNPrate$MinorFreq > 0.05]
selectedsnp <- sample(segregating, 2000)


geno.mat.sub <- geno.mat[selectedind, selectedsnp]
dim(geno.mat.sub)

metadata <- data.table(id    = indset,
                       SLiMid= scan(paste0(data.path, 
                                           "nonWF_st1_gen_6000_inds.txt")),
                       age   = scan(paste0(data.path, 
                                           "nonWF_st1_gen_6000_age.txt")),
                       strata= scan(paste0(data.path, 
                                           "nonWF_st1_gen_6000_sbpop.txt"), 
                                    what = "character"),
                       sex   = scan(paste0(data.path, 
                                           "nonWF_st1_gen_6000_sex.txt"), 
                                    what = "character"))
metadata.sub <- metadata[id %in% rownames(geno.mat.sub),]

# ---------------
# Perform the PCA
# ---------------

pca.slim     <- dudi.pca(geno.mat.sub, 
                         center=TRUE, scale=TRUE, nf = 30, scannf=FALSE)
pc.explained <- pca.slim$eig / sum(pca.slim$eig) * 100

# ---------------
# PCA plots
# ---------------


ggstruct <- pca.slim$li[, c(1,2,3)]
ggstruct <- merge(ggstruct, metadata.sub, by.x = "row.names", by.y = "id")

ggplot(ggstruct, aes(x=Axis1, y=Axis2, colour=factor(strata))) + 
  geom_point() + 
  theme_web_bw() + scale_colour_brewer(palette="Set2") +
  labs(x=paste0("PC1 (", signif(pc.explained[1], digits=3), "% of variance explained)"), 
       y=paste0("PC2 (", signif(pc.explained[3], digits=2), "% of variance explained)"), 
       color="Subpopulation")

ggplot(ggstruct, aes(x=Axis1, y=Axis3, colour=factor(strata))) + 
  geom_point() + 
  theme_web_bw() + scale_colour_brewer(palette="Set2") +
  labs(x=paste0("PC1 (", signif(pc.explained[1], digits=3), "% of variance explained)"), 
       y=paste0("PC3 (", signif(pc.explained[3], digits=2), "% of variance explained)"), 
       color="Subpopulation")
