# ==============================================================================
# Simple nearest-neighbour clustering function to assign points in kinship
# statistic space to a particular relatedness class.
# Author: Luke Lloyd-Jones 
# Date started: 24/03/2022
# Date updated: 24/03/2022
# ==============================================================================


clustKs <- function(ks, cent_sz = 3)
{
  # Inputs
  #   ks - the kinship statistics matrix from PCRelata
  #   cent_sz - size of centers 2 or 3?
  # Outouts
  #  clst - the pairs of individuals and their clustering
  #         from nearest neighbours
  kvals <- ks[, c("k0", "kin", "k2")]
  cents <- matrix(c(0,    0.25, 0.50,  0.75,   1, 
                    0.25, 0.25, 0.125, 0.063,  0,
                    0,    0.25, 0,     0,      0), 
                  nrow = 5, ncol = 3)
  cents <- cents[, 1:cent_sz]

  kd1 <- apply(kvals, 1, function(x) {sum((x - cents[1, ])^2)})
  kd2 <- apply(kvals, 1, function(x) {sum((x - cents[2, ])^2)})
  kd3 <- apply(kvals, 1, function(x) {sum((x - cents[3, ])^2)})
  kd4 <- apply(kvals, 1, function(x) {sum((x - cents[4, ])^2)})
  kd5 <- apply(kvals, 1, function(x) {sum((x - cents[5, ])^2)})
 
  kd   <- cbind(kd1, kd2, kd3, kd4, kd5)
  rel.set <- c("POP", "FSP", "2ndDeg", "3rdDeg", "U")
  clst <- apply(kd, 1, function(x) {rel.set[which.min(x)]})
  return(clst)
}