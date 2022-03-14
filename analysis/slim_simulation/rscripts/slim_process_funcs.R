# ==============================================================================
# Functions necessary to process SLiM output for looking at  
# Derived from PLODf_2al_tablefree_kira.R
# Authors: Kira Villiers
# Date started: 14/01/2021
# Date updated: 26/01/2021
# ==============================================================================

library(SNPRelate)

# Find kinships, as determined by an actual pedigree file, and optionally
# also find what PCRelate thresholding suggests kinships are. Produces
# a validation plot with actual kinships overlayed.

SLiM_validate <- function(PLOD, 
                          pcrelate = NULL, 
                          homefolder, 
                          longids = "nonWF_st1_gen_6000_inds.txt", 
                          subpops = "nonWF_st1_gen_6000_sbpop.txt", 
                          ages    = "nonWF_st1_gen_6000_age.txt",
                          ...) {
  # assumes pedigrees are in "mating_multi_gen_p1.txt" 
  # and  "mating_multi_gen_p2.txt"
  require(data.table)
  
  k <- mypcrel$kinBtwn
  
  # match index 'i' IDs to long 'id' IDs:
  ids <- data.table(longid=paste0("id", t(fread(paste0(homefolder, longids)))))
  ids$indexid <- paste0("i", 0:(length(ids$longid)-1))
  
  # ages
  ags <- data.table(ag=t(fread(paste0(homefolder, ages), header=FALSE)))
  ags$indexid <- paste0("i", 0:(length(ags$ag)-1))
  k$age1 <- ags$ag[match(k$ID1, ags$indexid)]
  k$age2 <- ags$ag[match(k$ID2, ags$indexid)]
  
  k$longID1 <- ids$longid[match(k$ID1, ids$indexid)]
  k$longID2 <- ids$longid[match(k$ID2, ids$indexid)]
  
  pop <- data.table(p=t(fread(paste0(homefolder, subpops), header=FALSE)))
  pop$indexid <- paste0("i", 0:(length(pop$p)-1))
  
  k$pop1 <- pop$p[match(k$ID1, pop$indexid)]
  k$pop2 <- pop$p[match(k$ID2, pop$indexid)]
  k$crosspop <- FALSE
  k$crosspop[k$pop1 != k$pop2] <- TRUE
 # k$pop1 <- NULL
 # k$pop2 <- NULL
  
  # The above 'i'd to long 'id' conversion is necessary because we 
  # only have the pedigree
  # in terms of long ids.
  ped.1 <- fread(paste0(homefolder,  "mating_multi_gen_p1.txt"))
  if (file.exists(paste0(homefolder, "mating_multi_gen_p2.txt"))) {
    ped.1 <- rbind(ped.1, fread(paste0(homefolder, "mating_multi_gen_p2.txt")))
  }
  colnames(ped.1) <- c("GenBorn", "mumid", "dadid", "id", "Sex")
  ped.1$mumid <- paste0("id", ped.1$mumid)
  ped.1$dadid <- paste0("id", ped.1$dadid)
  ped.1$id <- paste0("id", ped.1$id)
  
  # Work out actual kinships. For the moment it only 
  # works out POPs, FSPs and HSPs.
  # All others are labelled as unrelated.
  k2 <- findkinships(pedigree = ped.1, 
                     sampled = unique(c(k$longID1,k$longID2)))
  #ids$longid[match(kinferdata_g$info[,"Our_sample"],ids$indexid)]
  
  k <- rbind(merge(k, k2, by.x=c("longID1","longID2"), by.y=c("ID1","ID2")),
             merge(k, k2, by.x=c("longID1","longID2"), by.y=c("ID2","ID1")))
  
  
  if (!is.null(pcrelate)) {
    # Work out the kinships, as estimated by PCRelate:
    k3 <- k
    k3$status <- "U"
    k3$d1 <- k3$kin > 2^(-(1 + 3/2)) & k3$kin <= 2^(-(1 + 1/2)) #degree 1
    k3$d2 <- k3$kin > 2^(-(2 + 3/2)) & k3$kin <= 2^(-(2 + 1/2)) #degree 2
    k3$status[k3$d1 == TRUE] <- "FSP"
    k3$status[k3$d1 == TRUE & k3$k0 < 2^(-9/2)] <- "POP"
    # PCRelate cannot distinguish between HSP, aunt/niece, grandparent/grandchild. 
    # They are all in this category
    k3$status[k3$d2 == TRUE] <- "HSP" 
 
    
    k <- rbind(merge(k,  k3[,c("ID1","ID2","status")], 
                     by.x=c("ID1","ID2"), 
                     by.y=c("ID1","ID2")),
               merge(k,  k3[,c("ID1","ID2","status")], 
                     by.x=c("ID1","ID2"), 
                     by.y=c("ID2","ID1")))
    names(k)[names(k) == 'status'] <- 'estkinship.PC'
  }
  
  # attr(x = k, which="mean_UP")       <- attr(x = PLOD, which="mean_UP")
  # attr(x = k, which="mean_HSP")      <- attr(x = PLOD, which="mean_HSP")
  # attr(x = k, which="mean_FSP")      <- attr(x = PLOD, which="mean_FSP")
  # attr(x = k, which="mean_POP")      <- attr(x = PLOD, which="mean_POP")
  
  #plotPLOD(k, validations=TRUE, ...)
  return(k)
}

findkinships <- function(pedigree, sampled = NULL) {
  require(data.table)
  # pass me a pedigree file eg from SLiM -> for now, see the bottom of this script 
  # for an example. Will set out documentation soon. 
  
  # Check pedigree seems valid
  assertthat::assert_that(length(pedigree[["mumid"]]) > 0 && 
                          length(pedigree[["dadid"]]) > 0 && 
                          length(pedigree[["id"]])    > 0,
                          msg="`pedigree` must have columns `mumid`, `dadid`, and `id`.")
  
  # Create or check sampled seems valid
  if (is.null(sampled)) {
    # assume all individuals were sampled
    sampled <- unique(c(pedigree[["mumid"]], pedigree[["dadid"]], pedigree[["id"]]))
    
  } else if (!any(sampled %in% pedigree[["id"]])) {
    assertthat::assert_that(any(sampled %in% pedigree[["id"]]),
                            msg="No samples in `sampled` exist in `pedigree`.")
  } else if (!all(sampled %in% pedigree[["id"]])) {
    warning("At least some entries in `sampled` have no known pedigree in `pedigrees`.")
  }
  
  # Precalculate some booleans
  pedigree$me.insample  <- pedigree$id %in% sampled
  pedigree$mum.insample <- pedigree$mumid %in% sampled
  pedigree$dad.insample <- pedigree$dadid %in% sampled
  
  # make the result table
  ks <- data.table::data.table(cbind(t(combn(sampled, 2)), "U"))
  colnames(ks) <- c("ID1","ID2","kinship")
  
  # Start the kinship hunting:
  # Find Parent Offspring Pairs in our sample
  POP <- data.table::rbindlist(list(pedigree[, .(id,mumid)], pedigree[, .(id,dadid)]), 
                               use.names=FALSE)
  colnames(POP) <- c("V1", "V2")
  # The left join keeps only the ones in the sample
  ks[POP, on=.(ID1 == V1, ID2 == V2), kinship := "POP"]
  ks[POP, on=.(ID2 == V1, ID1 == V2), kinship := "POP"]
  #nPOP <- pedigree[me.insample & mum.insample & !dad.insample,.N] + 
  #  pedigree[me.insample & dad.insample & !mum.insample,.N] + 
  #  pedigree[me.insample & mum.insample & dad.insample,.N]
  nPOP <- sum(ks$kinship == "POP", na.rm=TRUE)
  
  
  # Find all Full Sibling Pairs
  FSP <- data.table::data.table(NULL)
  # all parental pairs that occur multiple times in pedigree table
  pgroups <- unique(pedigree[duplicated(pedigree, by=c("mumid","dadid")),
                                                     c("mumid","dadid")])
  
  for (ind in 1:dim(pgroups)[1]) {
    sibs <- pedigree[mumid == pgroups$mumid[ind] & dadid == pgroups$dadid[ind],]
    #if (dim(sibs)[1] > 1) { #at least two siblings exist, redundant check
    FSP <- rbind(FSP, t(combn(sibs$id, 2)))
    #}
  }
  colnames(FSP) <- c("V1", "V2")
  
  ks[FSP, on=.(ID1 == V1, ID2 == V2), kinship := "FSP"]
  ks[FSP, on=.(ID2 == V1, ID1 == V2), kinship := "FSP"]
  assertthat::assert_that(sum(ks$kinship == "POP", na.rm=TRUE) == nPOP, 
                          msg="Pedigree file seems wrong: found a pair that's 
                          apparently both parent/offspring and full siblings")
  nFSP <- sum(ks$kinship == "FSP", na.rm=TRUE)
  
  #Find half-sibling pairs. They have a parent and child in the sample but not both.
  HSP <- data.table(NULL)
  for (parent in pedigree[, sum(me.insample), by=.(mumid)][V1 > 1,mumid]) {
    sibs <- pedigree[mumid == parent & me.insample,]
    cmb <- combn(sibs$id,2)
    HSP <- rbind(HSP, t(cmb))
  }
  for (parent in pedigree[,sum(me.insample),by=.(dadid)][V1 > 1,dadid]) {
    sibs <- pedigree[dadid == parent & me.insample,]
    cmb <- combn(sibs$id,2)
    HSP <- rbind(HSP, t(cmb))
  }
  HSP <- HSP[!data.table(FSP), on=c("V1","V2")]
  ks[HSP, on=.(ID1 == V1, ID2 == V2), kinship := "HSP"]
  ks[HSP, on=.(ID2 == V1, ID1 == V2), kinship := "HSP"]
  nHSP <- sum(ks$kinship == "HSP", na.rm=TRUE)
  
  # --------------
  # Profiling GGPs
  # --------------
  dim(ks) # Cohort gap has to be 
  k.sub <- k[(k$kin > 0.1  & k$kin < 0.13) &
             (k$k0  > 0.45 &  k$k0 < 0.55), ]
  
  k2 <- rbind(merge(k, ks, by.x=c("longID1","longID2"), by.y=c("ID1","ID2")),
             merge(k, ks,  by.x=c("longID1","longID2"), by.y=c("ID2","ID1")))
  table(k2$kinship)
  
  k.sub <- k2[(k2$kin > 0.09  & k2$kin < 0.15) &
              (k2$k0  > 0.35 &  k2$k0 < 0.65)
               & (k2$kinship != "HSP") , ]
  
  max(k.sub$age1 - k.sub$age2)
  
  # id4538637 id4550946 i2230 i2682
  for (i in seq(1, dim(k.sub)[1]))
  {
    p1 <- pedigree[id == k.sub$longID1[i], ]
    p2 <- pedigree[id == k.sub$longID2[i], ]
  
    p1.mum <- pedigree[id == p1$mumid, ]
    p1.dad <- pedigree[id == p1$dadid, ]
  
    p1.ps <- c(p1$mumid, p1$dadid)
    p1.gs <- c(p1.mum$mumid, p1.mum$dadid, 
               p1.dad$mumid, p1.dad$dadid)
  
    p2.mum <- pedigree[id == p2$mumid, ]
    p2.dad <- pedigree[id == p2$dadid, ]
  
    p2.ps <- c(p2$mumid, p2$dadid)
    p2.gs <- c(p2.mum$mumid, p2.mum$dadid, 
               p2.dad$mumid, p2.dad$dadid)
    # Test
    if (sum(p1.gs %in% p2.ps) == 2 |
        sum(p2.gs %in% p1.ps) == 2)
    {
      print('Aunt/Neice relationship')
    }
  }
  
  
  return(ks)
}


# Do preprocessing (MAF filtering + LD filtering) of crocodile dataset.
#
# nsnps and maxsample are maximum sizes. if n.pcs is not provided 
# interactive use is enabled.
#
# returns a list. The first element is a filename where the filtered
# data is saved, the second a kinference dataframe, the third the PCRelate 
# output dataframe, and the last a ggplot object containing the plot of
# the first two PCs
SLiM_preprocess <- function(homefolder, maxsample=0, nsnps=0, n.pcs=0,
                            genodata="nonWF_st1_gen_6000.vcf", pca=TRUE) {
  require(SNPRelate)
  require(GWASTools)
  require(data.table)
  
  showfile.gds(closeall = TRUE)
  
  snpgdsVCF2GDS(vcf.fn=paste0(homefolder, genodata),
                out.fn=paste0(homefolder, "temp.gds"))
  
  showfile.gds(closeall = TRUE)
  
  f <- snpgdsOpen(paste0(homefolder, "temp.gds"))
  snpset   <- read.gdsn(index.gdsn(f, "snp.id"))
  indset   <- read.gdsn(index.gdsn(f, "sample.id"))
  
  #if (is.null(maxsample)) {
  if (maxsample <= 0) {
    selectedind <- indset
  } else {
    selectedind <- sample(indset, maxsample)
  }
  
  SNPrate <- snpgdsSNPRateFreq(f, with.snp.id = TRUE,
                               sample.id = selectedind)
  #hist(SNPrate$MinorFreq)
  segregating <- SNPrate$snp.id[SNPrate$MinorFreq > 0.05]
  
  ldset    <- snpgdsLDpruning(f, 
                              method    = "corr", 
                              ld.threshold = sqrt(0.1), 
                              verbose      = FALSE,
                              slide.max.n  = 500,
                              sample.id    = selectedind)
  pruneld <- unlist(ldset, use.names = FALSE)
  
  selectedsnp <- intersect(segregating, pruneld)
  
  cat("After preprocessing there are",  length(selectedind), "samples and",
      length(selectedsnp), "loci left in the dataset.\n")
  
  showfile.gds(closeall = TRUE)
  
  pcdat <- runPCRelate(homefolder, selectedind, selectedsnp, nsnps, n.pcs)
  mypcrelate  <- pcdat[[1]]
  selectedind <- pcdat[[2]]
  selectedsnp <- pcdat[[3]]
  snpnames.prefix <- pcdat[[4]] #you can really tell runPCRelate used to be 
  #part of this function but I last-minute lifted it out to be able to use 
  #it for the crocodile dataset too, sorry. Would make these two less 
  #weirdly intertwined if there was time
  
  showfile.gds(closeall = TRUE)
  
  mypcrelate$mu <- mypcrelate$mu[,paste0(snpnames.prefix,selectedsnp)]
  gdsSubset(parent.gds=paste0(homefolder, "temp.gds"),
            sub.gds=   paste0(homefolder, "tempsubset.gds"),
            sample.include=selectedind,
            snp.include=   selectedsnp
  )
  
  showfile.gds(closeall = TRUE)
  
  ggpca <- NULL
  if (pca == TRUE) {
    f2 <- snpgdsOpen(paste0(homefolder, "temp.gds"))
    geno.mat <- snpgdsGetGeno(f2, snpfirstdim = F)
    rownames(geno.mat) <- indset
    colnames(geno.mat) <- paste0(snpnames.prefix,snpset)
    
    geno.mat.sub <- geno.mat[selectedind,paste0(snpnames.prefix,selectedsnp)]
    
    metadata <- data.table(id    = indset,
                           SLiMid= scan(paste0(homefolder, "nonWF_st1_gen_6000_inds.txt")),
                           age   = scan(paste0(homefolder, "nonWF_st1_gen_6000_age.txt")),
                           strata= scan(paste0(homefolder, "nonWF_st1_gen_6000_sbpop.txt"), what="character"),
                           sex   = scan(paste0(homefolder, "nonWF_st1_gen_6000_sex.txt"), what="character"))
    metadata.sub <- metadata[id %in% rownames(geno.mat.sub),]
    
    #ggpca <- plotPCA(geno.mat.sub, metadata.sub)
  }
  showfile.gds(closeall = TRUE)
  
  return(list(paste0(homefolder, "tempsubset.gds"), mypcrelate))
}

# Called by the _preprocess functions. Calls the PCRelate functions
# from pcrelateModified_logistic.R
runPCRelate <- function(homefolder, selectedind, selectedsnp, nsnps, n.pcs) {
  # Then run through PCRelate processing
  require(BiocParallel)
  require(GWASTools)
  require(GENESIS)
  genofile <- snpgdsOpen(paste0(homefolder, "temp.gds"))
  snpnames <- read.gdsn(index.gdsn(genofile, "snp.id"))
  snpnames.prefix <- ifelse(!any(is.numeric(snpnames)), "", "m")
  
  ibd.robust <- snpgdsIBDKING(genofile, 
                              sample.id = selectedind, 
                              snp.id    = selectedsnp)
  #hist(ibd.robust$kinship, 1000)
  showfile.gds(closeall = TRUE)
  
  geno     <- GdsGenotypeReader(filename = paste0(homefolder, "temp.gds"))
  genoData <- GenotypeData(geno)
  
  KINGmat           <- ibd.robust$kinship
  colnames(KINGmat) <- as.character(ibd.robust$sample.id)
  rownames(KINGmat) <- as.character(ibd.robust$sample.id)
  mypcair <- GENESIS::pcair(genoData, 
                            kinobj = KINGmat, 
                            divobj = KINGmat,
                            sample.include  = selectedind,
                            snp.include = selectedsnp)
  if (n.pcs < 2) {
    plot(seq(1, length(mypcair$varprop)), mypcair$varprop)
    cat("(for non-interactive use, set parameter n.pcs)\n")
    uinput <- readline(prompt="Choose number of PCs to run PCRelate on (see plot): ")
    n.pcs <- as.integer(uinput)
    if (is.na(n.pcs) || n.pcs < 2) stop("Invalid choice of number of PCs.")
  }
  #plot(mypcair)
  
  # ======= Do the PCrelate ========
  showfile.gds(closeall = TRUE)
  geno     <- GdsGenotypeReader(filename = paste0(homefolder, "temp.gds"))
  genoData <- GenotypeData(geno)
  seq_genoData <- GenotypeBlockIterator(genoData,  snpInclude = selectedsnp)
  
  mypcrelate   <- pcrelate(seq_genoData, 
                           sample.include = selectedind,
                           ibd.probs      = TRUE,
                           maf.bound.method = "filter",
                           pcs            = mypcair$vectors[, 1:n.pcs], 
                           training.set   = mypcair$unrels,
                           small.samp.correct = TRUE)
  #colnames(mypcrelate$mu) <- paste0(snpnames.prefix,selectedsnp)
  cat("Finished calculating PCRelate values...\n")
  
  # filterset <- colSums(mypcrelate$mu, #>= 1 | mypcrelate$mu <= 0, 
  #                      na.rm = FALSE) #purpose: propagate NAs
  # filterbool <- !is.na(filterset) #!is.na(filterset) & filterset == 0
  # selectedsnp <- selectedsnp[filterbool]
  # 
  # if (nsnps > 0 & length(selectedsnp) > nsnps) {
  #   selectedsnp <- sample(selectedsnp, nsnps)
  # }
  # cat(length(selectedsnp), "SNPs have been kept")
  
  return(list(mypcrelate, selectedind, selectedsnp, snpnames.prefix))
}
