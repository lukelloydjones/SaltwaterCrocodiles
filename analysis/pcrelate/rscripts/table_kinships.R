# ==============================================================================
# Kin to comparisons tabling
# Author: Luke Lloyd-Jones
# Date started: 04/04/2022
# Date updated: 04/04/2022
# ==============================================================================
library(knitr)
library(Matrix)
kin <- read.csv(paste0("pcrelate/results_2/kinship_stats_pcrelate_", 
                       area, "plus_distance_all_models.csv"), 
                header = T)

head(kin)
table(kin$BIO1, kin$BIO2)

# Make the CMC coastal plains area the PCB

kin$BIO1[kin$BIO1 == "Coastal-Plains-CMC"] <- "Princess-Charlotte-Bay"
kin$BIO2[kin$BIO2 == "Coastal-Plains-CMC"] <- "Princess-Charlotte-Bay"

# Table again

tbl <- table(kin$BIO1, kin$BIO2)
topBot <- function(x)
{
  for (i in seq(1, dim(x)[1]))  
  {
  for (j in seq(i, dim(x)[2]))
  {
    if (i != j)
    {
      x[i, j] = as.numeric(x[i, j]) + as.numeric(x[j, i])
    }
  }
  }
  return(x)
}
tbl <- topBot(tbl)

tbl[lower.tri(tbl)] = t(tbl)[lower.tri(tbl)]

ord.nms <- c("Gulf-Plains-ALD", "Gulf-Plains-NFD", "Gulf-Plains-MGD",
             "North-West-Cape-York", "North-East-Cape-York",
             "Princess-Charlotte-Bay",   "Coastal-Plains-CA", 
             "Coastal-Plains-APrR", "Fitzroy")
tbl[lower.tri(tbl)] = t(tbl)[lower.tri(tbl)]
tbl <- tbl[ord.nms, ord.nms]
tbl[lower.tri(tbl, diag = F)] <- "--"

kable(tbl, format = 'latex')

# Make a template
tbl.tmp <- tbl
tbl.tmp[upper.tri(tbl.tmp, diag = T)] <- 0
tbl.tmp[lower.tri(tbl.tmp, diag = T)] <- 0
tbl.tmp <- apply(tbl.tmp, 2, as.numeric)
rownames(tbl.tmp) <- colnames(tbl.tmp)

# -----------
# Table POPs
# -----------

tbl.tmp.2 <- tbl.tmp
kin.pops <- kin[kin$ClstLSVMKinBig == "POP", ]
tbl.pops <- table(kin.pops$BIO1, kin.pops$BIO2)
tbl.tmp.2[rownames(tbl.pops), colnames(tbl.pops)] <- tbl.pops
tbl.pops <- tbl.tmp.2
tbl.pops <- topBot(tbl.pops)
tbl.pops[lower.tri(tbl.pops, diag = F)] <- "--"


# -----------
# Table FSPs
# -----------

tbl.tmp.2 <- tbl.tmp
kin.fsps <- kin[kin$ClstLSVMKinBig == "FSP", ]
tbl.fsps <- table(kin.fsps$BIO1, kin.fsps$BIO2)
tbl.tmp.2[rownames(tbl.fsps), colnames(tbl.fsps)] <- tbl.fsps
tbl.fsps <- tbl.tmp.2
tbl.fsps <- topBot(tbl.fsps)
tbl.fsps[lower.tri(tbl.fsps, diag = F)] <- "--"


# -----------
# Table 2ND
# -----------

tbl.tmp.2 <- tbl.tmp
kin.2nd <- kin[kin$ClstLSVMKinBig == "2nd-Deg", ]
tbl.2nd <- table(kin.2nd$BIO1, kin.2nd$BIO2)
tbl.tmp.2[rownames(tbl.2nd), colnames(tbl.2nd)] <- tbl.2nd
tbl.2nd <- tbl.tmp.2
tbl.2nd <- topBot(tbl.2nd)
tbl.2nd[lower.tri(tbl.2nd, diag = F)] <- "--"

tbl.kin <- tbl.2nd
tbl.kin[] <- paste(tbl.pops, tbl.fsps, tbl.2nd, sep=",")
tbl.kin[lower.tri(tbl.kin, diag = F)] <- "--"
tbl.kin


# Try and table the POPS

mta.fil.pth <- "croc_dart_round_2/meta_data_derived/croc_info_mapped_bioregion_updated_3.csv"
mta.dta <- read.csv(mta.fil.pth, header = T)
rownames(mta.dta) <- paste0("IND", mta.dta$Sample_Id)

mta.dta.I1 <- mta.dta[kin.pops$IDS1, ]
mta.dta.I1 <- mta.dta.I1[, c("Sample_Id", "Length_Class", 
               "Sex_MF", "Capture_Location", "Lat", "Long", "Date_Sampled_2",
               "Bioregion_New_Nm")]
colnames(mta.dta.I1) <- paste0("IND1_", colnames(mta.dta.I1))

mta.dta.I2 <- mta.dta[kin.pops$IDS2, ]
mta.dta.I2 <- mta.dta.I2[, c("Sample_Id", "Length_Class", 
                    "Sex_MF", "Capture_Location", "Lat", "Long", "Date_Sampled_2",
               "Bioregion_New_Nm")]
colnames(mta.dta.I2) <- paste0("IND2_", colnames(mta.dta.I2))

pop.etra <- cbind(kin.pops, mta.dta.I1, mta.dta.I2)

pop.etra <- pop.etra[,c("IDS1", "IDS2", "kin", "k0", "k2",                     
                        "ClstKin", "ClstLSVMKinBig",                            
                        "BIO1",  "BIO2",                            
                        "IND1_Sex_MF" , "IND2_Sex_MF",                
                        "IND1_Capture_Location", "IND2_Capture_Location", "DIST_1T2",  
                        "IND1_Date_Sampled_2",   "IND2_Date_Sampled_2",      
                        "IND1_Length_Class", "IND2_Length_Class")]
pop.etra$DAYS_DIFF <- abs(as.Date(pop.etra$IND1_Date_Sampled_2) - 
                       as.Date(pop.etra$IND2_Date_Sampled_2))

write.csv(pop.etra, "pcrelate/results_2/pops_classified.csv",
          row.names = F)

