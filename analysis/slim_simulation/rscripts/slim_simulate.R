# ==============================================================================
# SCript to generate the simulation results for getting a handle on false
# positive rates and false negative rates from PCRelate
# Authors: Luke Lloyd-Jones and Kira Villiers
# Date started: 07/03/2022
# Date updated: 07/03/2022
# ==============================================================================
library(ggplot2)
library(ggwebthemes)
library(RColorBrewer)

# Class close
showfile.gds(closeall = TRUE)

condition <- "0pt0001pc"
srcfolder <- paste0("~/Desktop/slim/slim_sim_2/slim_mg_2p_migr_", condition, "/")
out.fl    <- "slim_simulation/out/"

prep <- SLiM_preprocess(srcfolder, maxsample = 500, nsnps = 3000, n.pcs = 3)

mypcrel <- prep[[2]]

k <- SLiM_validate(mypcrel, homefolder = srcfolder, pcrelate = T)

saveRDS(k,   paste0(out.fl, "pcrelate_results", condition, ".rds"))

colrs <- c("black", "orange", "gold", "blue")
names(colrs) <- c("U","HSP","FSP","POP")

g <- ggplot(k, 
            aes(x = k0, y = kin, colour = factor(kinship))) +
            geom_point() + 
            theme_web_bw() + 
            scale_colour_manual(values = colrs) +
            geom_hline(yintercept = 2^(-(1 + 1/2)), 
                       colour = colrs["U"],   
                       linetype="dashed") +
            geom_hline(yintercept = 2^(-(2 + 1/2)), 
                       colour = colrs["POP"], 
                       linetype="dashed") +
            geom_hline(yintercept = 2^(-(2 + 3/2)), 
                       colour = colrs["HSP"], 
                       linetype="dashed") +
            theme(legend.position = c(0.85, 0.8),
                  legend.title = element_text(size = 20, face = "bold"),
                  legend.text  = element_text(size = 20, face = "bold"),
                  axis.text    = element_text(size = 20, face = "bold"),
                  axis.title.x = element_text(size = 22, face = "bold"),
                  axis.title.y = element_text(size = 22, face = "bold")) +
            guides(colour=guide_legend(title = "Kinship")) +
            xlab("k(0)") + 
            ylab("Kinship coefficient")

g.1pc   <- g
g.01pc  <- g
g.001pc <- g

k.hsp <- k[k$kinship == "HSP", ]
mean(k.hsp$kin)
mean(k.hsp$k0)

k.pop <- k[k$kinship == "POP", ]
mean(k.pop$kin)





