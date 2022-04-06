# ==============================================================================
# PCrelate for crocodiles
# Author: Luke Lloyd-Jones
# Date started: 08/03/2022
# Date updated: 29/03/2022 
# ==============================================================================
library(GENESIS)
library(SNPRelate)
library(dplyr)
library(GWASTools)
library(adegenet)
library(hierfstat)
library(ggplot2)
library(ggwebthemes)
library(gridExtra)
library(geosphere)
library(dartR)
library(e1071)

# Read croc info

showfile.gds(closeall = TRUE)
in.pth <- "radiator/kira_radiator_filtered/crocs_qced_genid_allalles.Rdata"

raddata  <- get(load(in.pth))
crocs.gl <- gi2gl(raddata, parallel = FALSE, verbose = T)

# Use the DARTR functions to write out

#gl2gds(crocs.gl, 
#       outpath = ".",
#       outfile = "crocs_qced_genid_allalles.gds")

showfile.gds(closeall = TRUE) # Here if you need

# Read the GDS now
gds.pth <- "radiator/kira_radiator_filtered/crocs_qced_genid_allalles.gds"
genofile <- snpgdsOpen(gds.pth)

# IDs and bioregions

samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
srta    <- raddata@strata 
all(samp.id == srta$INDIVIDUALS)

# Subset to groups of bioregions
area <- 'west_coast'
if (area == 'west_coast')
{
  srta.west <- srta[(srta$STRATA == "2__North-West-Cape-York" | 
                     srta$STRATA == "1c__Gulf-Plains-NFD"   | 
                     srta$STRATA == "1b__Gulf-Plains-ALD"   | 
                     srta$STRATA == "1d__Gulf-Plains-MGD"), ]
  srta <- srta.west
} 

area <- "all"

# Prune markers by linkage disequilibrium

snpset   <- snpgdsLDpruning(genofile,
                            method       = "corr", 
                            ld.threshold = sqrt(0.1), 
                            verbose      = FALSE,
                            slide.max.n  = 500, 
                            sample.id    = srta$INDIVIDUALS)

pruned <- unlist(snpset, use.names = FALSE)
length(pruned)


# Subset on missing rate and allele frequency
RV        <- snpgdsSNPRateFreq(genofile, 
                               snp.id = pruned, 
                               sample.id   = srta$INDIVIDUALS,
                               with.snp.id = TRUE)
af.df      <- data.frame(RV)
af.df.pt05 <- af.df[RV$MinorFreq > 0.05 , ]

# ---------------------------
# Begin the PCRelate pipeline
# ---------------------------
# Compute the IBD robust statistics 

ibd.robust <- snpgdsIBDKING(genofile, 
                            sample.id = srta$INDIVIDUALS,
                            snp.id = af.df.pt05$snp.id)

plt <- hist(ibd.robust$kinship, 35000, plot = F)
lg <- log(plt$counts)
lg[is.infinite(lg)] <- NA
jpeg(paste0("pcrelate/results_2/king_robust_crocs_", area, ".jpeg"),
     pointsize = 26, 
     quality = 300, 
     bg = "white", res = NA, width = 1300, height = 800)
plot(plt$breaks[-1], lg, 
     type = 'h', 
     lwd  = 10, 
     lend = 2, 
     xlim = c(-1.5, 0.5), 
     xlab = "KING-robust kinship statistic",
     ylab = "log(counts)")
abline(v = 0.125, col = 'red')
abline(v = 0.25,  col = 'red')
dev.off()


rowMeans(ibd.robust$kinship)

showfile.gds(closeall = TRUE)

# ------------------------------
# Move on to the PCAir component
# ------------------------------

geno     <- GdsGenotypeReader(filename = "radiator/kira_radiator_filtered/crocs_qced_genid_allalles.gds")
genoData <- GenotypeData(geno)


KINGmat <- ibd.robust$kinship
colnames(KINGmat) <- as.character(ibd.robust$sample.id)
rownames(KINGmat) <- as.character(ibd.robust$sample.id)
mypcair <- pcair(genoData, 
                 kinobj = KINGmat, 
                 divobj = KINGmat,
                 sample.include = srta$INDIVIDUALS,
                 snp.include = af.df.pt05$snp.id)
summary(mypcair)


jpeg(paste0('pcrelate/results_2/pcair_pc1_v_pc2_prop_var_', area ,'.jpeg'),
     pointsize = 26, 
     quality = 300, 
     bg = "white", res = NA, width = 2300, height = 1300)
par(mfrow = c(1, 2))
plot(mypcair)
plot(seq(1, length(mypcair$varprop)), mypcair$varprop, 
     xlab = "PC number", 
     ylab = "Proportion variance explained by PC",
     xlim = c(0, 50))
dev.off()

# ------------------------------
# Back to the PCRelate statistic
# ------------------------------

seq_genoData <- GenotypeBlockIterator(genoData, 
                                      snpInclude = af.df.pt05$snp.id)
mypcrelate   <- pcrelate(seq_genoData, 
                         sample.include     = srta$INDIVIDUALS,
                         pcs                = mypcair$vectors[, 1:10], 
                         training.set       = mypcair$unrels,
                         ibd.probs          = TRUE,
                         maf.bound.method   = "filter",
                         small.samp.correct = TRUE)

# ------------------------------------------
# Order the results by kinship statistics
# ------------------------------------------

kinpairs   <- mypcrelate$kinBtwn
kinpairs.o <- kinpairs[order(kinpairs$kin, decreasing = T), ]

write.csv(kinpairs.o, 
          paste0("pcrelate/results_2/kinship_stats_pcrelate_", area, ".csv"),
          row.names = F, 
          quote = T)

# --------------------------------------
# Let's screen for unrelated individuals
#  - this was done using a kin.thresh - 0.05
# --------------------------------------
srta.unrel <- srta[mypcair$unrels, ]
write.csv(srta.unrel, 
          "samples_edna/unrelated_set.csv",
          row.names = F, quote = T)

# -------------------------------------------------
# Log counts plot on kinship coefficient similar to
# Kinference
# -------------------------------------------------
kinpairs.o <- read.csv(paste0("pcrelate/results_2/kinship_stats_pcrelate_", 
                              area, ".csv"), 
                       header = T)

plt <- hist(kinpairs.o$kin, 35000, plot = F)
length(mypcrelate$kinBtwn$kin)

lg <- log(plt$counts)
lg[is.infinite(lg)] <- NA
jpeg(paste0("pcrelate/results_2/kinship_stat_hist_crocs_", area, ".jpeg"),
     pointsize = 26, 
     quality = 300, 
     bg = "white", res = NA, width = 1300, height = 800)
plot(plt$breaks[-1], lg, 
     type = 'h', 
     lwd  = 10, 
     lend = 2, 
     xlim = c(-0.1, 0.4), 
     xlab = "PC relate kinship statistic",
     ylab = "log(counts)")
abline(v = 0.125, col = 'red')
abline(v = 0.25,  col = 'red')
dev.off()

# -------------------------------------------------
# Kinship statistics versus k0 
# -------------------------------------------------
kinpairs.o <- read.csv(paste0("pcrelate/results_2/kinship_stats_pcrelate_allplus_distance.csv"), 
                       header = T)

# Do the clustering

# Euclid
# ------
kinpairs.o$ClstKin     <- clustKs(kinpairs.o)

# SVM
# ---
dat <- data.frame(kinpairs.o[, c('kin', 'k0', 'k2')])

# Diff rates

s     <- "0pt001pc_litter"
svm.s <- get(load(paste0("slim_simulation/svm_stacks/svm_linear", s, ".Rda")))
new <- predict(svm.s, dat)
kinpairs.o$ClstLSVMKin001 <- new

s <- "0pt01pc_litter"
svm.s <- get(load(paste0("slim_simulation/svm_stacks/svm_linear", s, ".Rda")))
new <- predict(svm.s, dat)
kinpairs.o$ClstLSVMKin01 <- new

s <- "0pt0001pc_litter"
svm.s <- get(load(paste0("slim_simulation/svm_stacks/svm_linear", s, ".Rda")))
new <- predict(svm.s, dat)
kinpairs.o$ClstLSVMKin0001 <- new

# All together
svm.s <- get(load("slim_simulation/svm_stacks/svm_linear_all_litter.Rda"))
new <- predict(svm.s, dat)
kinpairs.o$ClstLSVMKinBig <- new

table(kinpairs.o$ClstLSVMKin01)
table(kinpairs.o$ClstLSVMKin001)
table(kinpairs.o$ClstLSVMKin0001)
table(kinpairs.o$ClstLSVMKinBig)

kin.pops <- kinpairs.o[kinpairs.o$ClstLSVMKinBig == "POP", ]
mean(kin.pops$kin);mean(kin.pops$k0);mean(kin.pops$k2);
max(kin.pops$kin); max(kin.pops$k0); max(kin.pops$k2)

kin.fsps <- kinpairs.o[kinpairs.o$ClstLSVMKinBig == "FSP", ]
round(c(mean(kin.fsps$kin), mean(kin.fsps$k0), mean(kin.fsps$k2)), 3)
round(c(min(kin.fsps$kin),  min(kin.fsps$k0),  min(kin.fsps$k2)), 3)

kin.hsps <- kinpairs.o[kinpairs.o$ClstLSVMKinBig == "2nd-Deg", ]
round(c(mean(kin.hsps$kin), mean(kin.hsps$k0), mean(kin.hsps$k2)), 3)
round(c(min(kin.hsps$kin),  min(kin.hsps$k0),  min(kin.hsps$k2)), 3)

# ----
# Draw
# ----

# Source the plotting functions
source("slim_simulation/rscripts/plot_funcs.R")
source("slim_simulation/rscripts/cluster_func.R")

# Just a subsamp
r.smp <- sample(seq(1, dim(kinpairs.o)[1]), 25000)
kinpairs.o.2 <- kinpairs.o[r.smp, ]

kinpairs.o.2 <- kinpairs.o
kinpairs.o.2$kinship <- kinpairs.o.2$ClstKin

kin.g.t  <- drawkink0(kinpairs.o.2, clstr = 1)
k0k2.g.t <- drawk2k0(kinpairs.o.2,  clstr = 1)

kinpairs.o.2$kinship <- kinpairs.o.2$ClstLSVMKinBig

kinpairs.o.2$kinship <- droplevels(kinpairs.o.2$kinship)
levels(kinpairs.o.2$kinship) <- c("FSP", "2nd-Deg", "POP", "U")

kin.g.c  <- drawkink0(kinpairs.o.2, clstr = 3)
k0k2.g.c <- drawk2k0(kinpairs.o.2,  clstr = 3)

png(filename = paste0("pcrelate/results_2/all_clustering.png"),
    bg = "white", width =  18, height = 16, 
    units = 'in', res = 300)
grid.arrange(kin.g.t, k0k2.g.t,
             kin.g.c, k0k2.g.c, ncol = 2)
dev.off()


# ---------------------------------------------
# Distance plot
# ---------------------------------------------

mta.fil.pth <- "croc_dart_round_2/meta_data_derived/croc_info_mapped_bioregion_updated_3.csv"
mta.dta <- read.csv(mta.fil.pth, header = T)

# Match these with the unrelated set
rownames(mta.dta) <- paste0("IND", mta.dta$Sample_Id)
mta.dta.lt.lg     <- mta.dta[, c("Long", "Lat")]

kinpairs.o$IDS1 <- paste0("IND", gsub("-.*", "", kinpairs.o$ID1))
kinpairs.o$IDS2 <- paste0("IND", gsub("-.*", "", kinpairs.o$ID2))

lt.lg.1           <- mta.dta.lt.lg[kinpairs.o$IDS1, ]
colnames(lt.lg.1) <- paste0(c("Lat", "Long"), "_ID1")
lt.lg.2           <- mta.dta.lt.lg[kinpairs.o$IDS2, ]
colnames(lt.lg.2) <- paste0(c("Lat", "Long"), "_ID2")

lt.lg.1.m <- as.matrix(lt.lg.1)
lt.lg.2.m <- as.matrix(lt.lg.2)

dst.btw <- data.frame(DIST = rep(NA, dim(lt.lg.2.m)[1]))
for (i in seq(1, dim(dst.btw)[1]))
{ 
  if (i %% 100 == 0) {print(i)}
  dst.btw[i, ] <- distm(lt.lg.1.m[i, ], 
                        lt.lg.2.m[i, ], 
                        fun = distGeo)
}

# Join

kinpairs.o$DIST_1T2 <- dst.btw$DIST / 1000

mta.dta.bio <- mta.dta[, c("Bioregion_New", "Bioregion_New_Nm")]
kinpairs.o$BIO1 <- mta.dta.bio[kinpairs.o$IDS1, "Bioregion_New_Nm"]
kinpairs.o$BIO2 <- mta.dta.bio[kinpairs.o$IDS2, "Bioregion_New_Nm"]

area <- "all"
write.csv(kinpairs.o, 
          paste0("pcrelate/results_2/kinship_stats_pcrelate_", 
                 area, "plus_distance_all_models.csv"),
          row.names = F, 
          quote = T)

# ------------------- 
# Distance plot
# -------------------

colrs        <- c("black", "orange", "gold", "blue")
names(colrs) <- c("U", "2nd-Deg", "FSP", "POP")

# Oragnise levels 
kinpairs.o.smp <- kinpairs.o[sample(seq(1, dim(kinpairs.o)[1]), 10000), ]
kinpairs.o.smp$ClstLSVMKin <- droplevels(kinpairs.o.smp$ClstLSVMKin)
levels(kinpairs.o.smp$ClstLSVMKin) <- c("FSP", "2nd-Deg", "POP", "U")


kinpairs.o$ClstLSVMKinBig <- droplevels(kinpairs.o$ClstLSVMKinBig)
levels(kinpairs.o$ClstLSVMKinBig) <- c("FSP", "2nd-Deg", "POP", "U")
kinpairs.o$ClstLSVMKin <- kinpairs.o$ClstLSVMKinBig

gdist <- drawDistKin(kinpairs.o)
gdist

png(filename = paste0("pcrelate/results_2/clustering_plus_dist.png"),
    bg = "white", width =  18, height = 16, 
    units = 'in', res = 300)
 grid.arrange(kin.g.c, k0k2.g.c,
             gdist,
              layout_matrix = rbind(c(1, 2), c(3, 3)),
              ncol = 2)
dev.off()


# -------------------
# Distance summaries
# -------------------

kin.ps <- kinpairs.o[kinpairs.o$ClstLSVMKinBig != "U", ]
dim(kin.ps)

# Which had big distances?

kin.ps.bg.d <- kin.ps[kin.ps$DIST_1T2 > 800, ]
table(kin.ps.bg.d$BIO1, kin.ps.bg.d$BIO2)

# Pops at big distances

kin.ps.bg.d.pop <- kin.ps[kin.ps$DIST_1T2 > 400 &
                          kin.ps$ClstLSVMKinBig == "POP", ]

# Proportion of kin less than 50 km
sum(kin.ps$DIST_1T2 < 50) / length(kin.ps$DIST_1T2)

























