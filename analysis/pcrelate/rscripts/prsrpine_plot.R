# ==============================================================================
# Can we get a nice plot of the proserpine points and relationships?
# Author: Luke Lloyd-Jones
# Date started: 01/04/2022
# Date updated: 01/04/2022 
# ==============================================================================
library(ggmap)
library(ozmaps)
library(maps)
library(sp)
library(raster)
library(ggrepel)
library(RColorBrewer)
library(grid)
library(maptools)
# Propserpine set

mta.fil.pth <- "croc_dart_round_2/meta_data_derived/croc_info_mapped_bioregion_updated_3.csv"
mta.dta <- read.csv(mta.fil.pth, header = T)
rownames(mta.dta) <- paste0("IND", mta.dta$Sample_Id)


table(mta.dta$Bioregion_New_Nm)
mta.dta.prsp <- mta.dta[mta.dta$Bioregion_New_Nm == "Coastal-Plains-APrR", ]
mta.dta.prsp <- mta.dta.prsp[mta.dta.prsp$Place != "Sarina", ]

# La points

pts <- mta.dta.prsp

# La relativi

kin <- read.csv(paste0("pcrelate/results_2/kinship_stats_pcrelate_", 
                        area, "plus_distance.csv"), header = T)
kin$Long1 <- mta.dta[kin$IDS1, "Long"]
kin$Lat1  <- mta.dta[kin$IDS1, "Lat"] 
kin$Long2 <- mta.dta[kin$IDS2, "Long"]
kin$Lat2  <- mta.dta[kin$IDS2, "Lat"] 

kin.pros <- kin[(kin$BIO1 == "Coastal-Plains-APrR" & 
                 kin$BIO2 == "Coastal-Plains-APrR" & 
                 kin$ClstLSVMKin != "U"), ]
table(kin.pros$ClstLSVMKin)

kin.pros <- kin.pros[1:10, ]
kin.set  <- kin.pros

pp <- plotIndTlmtry(pts = pts, kin.set = kin.pros)

png(filename = 'pcrelate/results_2/proserpine_spider.jpeg',
          bg = "white", 
          width =  18, 
          height = 16, 
          units = 'in', 
         res = 300)
  print(pp)
dev.off()

plotIndTlmtry <- function(pts            = NULL,
                          kin.set        = NULL,                          
                          pth.t.ky       = "pcrelate/data/google_maps_api.txt", 
                          cty.sz         = 8,
                          pop.thr        = 3000,
                          lgs.pth        = "pcrelate/data/",
                          cpt            = NULL,
                          ggl.scl        = 1,
                          offsetmin      = 0.1,
                          offsetmax      = 0.1,
                          str.end.lbl.sz = 5,
                          plt.dta.pth    = "plot_data/",
                          plt.wtr.plys   = "pcrelate/data/water_polygons_big.Rda",
                          rvrs.pth       = "pcrelate/data/rivers_big.Rda",
                          rvrs.plt.sz    = 4,
                          rvrs.plt.thr   = 4,
                          bndry.dta.pth  = "pcrelate/data/",
                          tit.plt        = "Proserpine Region",
                          lg.pos         = c(0.93, 0.17))
{
  
  # Arguments in
  #                Default
  #  pts     = NULL,  a data frame containing telemetry points. Has to 
  #                       have 'Latitude' and 'Longitude' as two of the column 
  #                       names
  #  pth.t.ky    = NULL,  path to Google API key
  #  cty.sz      = 9,     a scalar to add to log values to control the size of
  #                       the city name plotted.
  #  pop.thr     = 35000, the population size in which to threshold the 
  #                       inclusion of certain cities.
  #  lgs.pth     = NULL,  path to logo files.
  #  cpt         = NULL,  an array of 2 elements with the center point in 
  #                       (long, lat) which to pull the Google maps plots. The
  #                       Google maps plot extracted with be centered at 'cpt'.
  #  ggl.scl     = 2,     this controls the difference in the scale that you
  #                       pull the google maps image. 
  #  offsetmin   = 1,     extend the bounding box of the telemetry data set to
  #                       make sure all the points fit in and there is enough
  #                       space to fit the logos and other legends etc.
  #  offsetmax   = 1,     Same as above. The units are lat long degrees.
  #  mid.nts.lbl.num = 5, 
  #  mid.nts.lbl.sz  = 3,
  #  arr.no.ctrl,
  #  water.not       = T,  USEFUL! Plotting all the data can be slow and if 
  #                       you are tuning parameters is just too slow. 
  #                       Set it to 5,000 if you want to plot things quickly
  #                       and tune them.
  #  plt.dta.pth = NULL,  Path to the Australia and state boundary shape file
  #  scl.var     = 10,    Some control over the size of the North arrow and 
  #                       scale bar
  #  lg.var      = 15,    Some control over the scale of the logos
  #  lg.pos      = c(0.93, 0.17) Vector of x, y coords positioning the SP legend   
  #  lg.cs.xmn   = 3      Deviation from xmin to start the CSIRO logo
  #  lg.cs.xmx   = 6      Deviation from xmin to end the CSIRO logo
  #  lg.cs.ymx   = 1      Deviation from ymin for height of CSIRO logo
  #  lg.cw.xmn   = 6.5    Deviation from xmin to end the CEWO logo
  #  lg.cw.xmx   = 12     Deviation from xmin to end the CEWO logo
  #  lg.cs.ymx   = 1      Deviation from ymin for height of CEWO logo
  #  scl.br.x.bf = 5      Buffer in degrees from from xmax to draw scale bar
  #  scl.br.y.df = 0.5    Buffer in degrees from from ymin to draw scale bar
  #  scl.dLon    = 100    Scale bar distance longitude 
  #  scl.dLat    = 54     Scale bar distance latitude 
  #  scl.dLeg    = 70     Scale bar distance length
  #  scl.dLgsz   = 7      Scale bar legend size
  #  scl.arwLgt  = 90     Scale bar arrow length
  #  scl.arwDst  = 114    Scale bar arrow distance
  #  scl.arwNsz  = 11     Scale bar north size
  # 
  #  Return
  #  pp - a ggplot object with the plot of all telemetry paths for the points
  #       passed to the function.
  
  # God-awful GOOGLE-API key. You needs this, which is a bit cheeky of google
  # I think you get a lot of 'free' google maps queries, which requires you
  # to register a credit card (Yuk). See here to do that
  # https://developers.google.com/maps/documentation/embed/get-api-key
  # Setup billing notifications and limits as well just in case.
  
  key_gmp = as.character(read.table(pth.t.ky))
  register_google(key = key_gmp)
  
  # All the cities of Australia
  aus.cty <- world.cities[world.cities$country.etc == "Australia", ]
  
  pts.sp <- SpatialPointsDataFrame(pts[, c("Long", "Lat")], pts)
  
  cpt <- coordinates(as(extent(pts.sp), "SpatialPolygons")) 
  
  # ------------------------
  # Basic map
  # ------------------------
  
  pts.loc <- extent(pts.sp)
  box     <- make_bbox(Longitude, 
                       Latitude, 
                       data = pts, 
                       f = 0.2)
  
  gmap <- get_googlemap(c(cpt),
                        zoom = calc_zoom(box) - ggl.scl,
                        scale = 2,
                        maptype = 'satellite')
  p  <- ggmap(gmap)
  
  box2 <- pts.loc
  box2[c(1, 3)] <- box2[c(1, 3)] - offsetmin
  box2[c(2, 4)] <- box2[c(2, 4)] + offsetmax
  
  # ----------------------------------------
  # Start, end and five midnights in between
  # ----------------------------------------
  
  pts$Label     <- NA
  pts$LabelSize <- pts$Length_Class / 3

  # -----------------------------------------------
  # Make the colours for different times of the day
  # -----------------------------------------------
  
  pts$Cols <- "green"
  
  # ---------------------------
  # Load up your rivers
  # ---------------------------
  
  shp.all.bigs <- get(load("shapefiles/all_big_rivers.shp"))
  dat_SP_LL <- spTransform(shp.all.bigs, 
                           CRS("+proj=longlat +datum=WGS84"))
  shp.all.bigs.frt <- fortify(dat_SP_LL)

  
  # Get city information
  
  aus.cty$SZ_CTY <- log(aus.cty$pop / max(aus.cty$pop)) + cty.sz
  aus.cty.35k    <- aus.cty[aus.cty$pop > pop.thr, ]
  
  # -------------
  # Get the logos
  # -------------
  
  lg.csiro <- get_png(paste0(lgs.pth, "CSIRO_Logo.svg.png"))
  #lg.cew   <- get_png(paste0(lgs.pth, "CEWO_logoStacked_HighRes.png"))
  
  # -----------------------
  # Adding the inset
  # -----------------------
  shape.sub <- get(load(paste0(bndry.dta.pth, "aus_coastline.Rda")))
  
  temp  <- data.frame(long = c(box2[1], box2[2], 
                               box2[2], box2[1], box2[1]),
                      lat  = c(box2[3], box2[3], box2[4], box2[4], 
                               box2[3]))
  g2 <- ggplotGrob(
    ggplot(data = shape.sub) +  
      geom_line(aes(x = long, y = lat, group = group), 
                fill = NA, colour = "black") + 
      theme_web_bw() +
      geom_path(data = temp, 
                aes(x = long, y = lat), 
                size = 1, colour = "red") +
      theme(axis.title.x = element_blank(),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y  = element_blank(),
            axis.ticks.y = element_blank())
  )
  
  # --------------------------------------
  # Compute the logo and scale positioning
  # --------------------------------------
  
  xdst <- box2[2] - box2[1]
  ydst <- box2[4] - box2[3]
  
  ist.cs.xmx <- xdst * 0.12
  ist.cs.ymx <- ydst * 0.2
  
  lg.cs.xmn <- xdst * 0.06
  lg.cs.xmx <- xdst * 0.1
  lg.cs.ymx <- ydst * 0.05
  
  lg.cw.xmn <- xdst * 0.18
  lg.cw.xmx <- xdst * 0.37
  lg.cw.ymx <- ydst * 0.05
  
  scl.br.x.bf <- xdst * 0.15
  scl.br.y.bf <- xdst * 0.03
  
  scl.dLon <- floor(xdst * 0.08 * 11.1)  * 10
  if (scl.dLon == 0)
  {
    scl.dLon <- 3
  }
  scl.dLat <- floor(ydst * 0.025 * 11.1) * 10
  if (scl.dLat == 0)
  {
    scl.dLat <- 2
  }
  scl.dLeg <- scl.dLat * 1.6
  scl.dLgsz  <- 5
  
  scl.arwLgt <- scl.dLat * 1.9
  scl.arwDst <- scl.dLat * 2.2
  scl.arwNsz <- 9
  
  ggl.tg.sz <- 4
  
  # -----------------------
  # Divine the plot
  # -----------------------
  adjustcolor("orange", alpha.f = 0.2)
  colrs        <- c("black", 
                    adjustcolor("orange", alpha.f = 0.1), 
                    adjustcolor("yellow", alpha.f = 0.4), 
                    adjustcolor("blue",   alpha.f = 0.5))
  names(colrs) <- c("U", "HSP", "FSP", "POP")
  title <- tit.plt

  
  pp <- p + 
    scale_colour_manual(values = colrs) +
    scale_x_continuous(limits = c(box2[1], box2[2]), expand = c(0, 0)) +
    scale_y_continuous(limits = c(box2[3], box2[4]), expand = c(0, 0)) + 
    geom_line(data = shp.all.bigs.frt, size = 1, 
              mapping = aes(x = long, y = lat, group = group), 
              colour  = brewer.pal(7, "Blues")[4]) +
    geom_segment(aes(x    = Long1, 
                     y    = Lat1, 
                     xend = Long2, 
                     yend = Lat2, color = ClstLSVMKin), 
                 data = kin.set) +
    geom_point(aes(x = Long, 
                  y = Lat), 
               data   = pts, 
               size   = pts$LabelSize, 
               colour = pts$Cols) +
    theme_web_bw() + 
    theme(legend.position = "none",
          legend.title = element_text(size = 20, face = "bold"),
          legend.text  = element_text(size = 20, face = "bold"),
          axis.text    = element_text(size = 20, face = "bold"),
          axis.title.x = element_text(size = 22, face = "bold"),
          axis.title.y = element_text(size = 22, face = "bold")) +
    guides(colour=guide_legend(title = "Kinship")) +
    xlab("Longitude") + 
    ylab("Latitude")  +
    inset(grob = g2, 
          xmin = box2[2] - ist.cs.xmx, 
          xmax = box2[2], 
          ymin = box2[4] - ist.cs.ymx, 
          ymax = box2[4]) +
    inset(lg.csiro, 
          xmin = box2[1] + lg.cs.xmn,
          xmax = box2[1] + lg.cs.xmx,  
          ymin = box2[3], 
          ymax = box2[3] + lg.cs.ymx)  + 
    geom_label(y = -Inf, x =  -Inf, 
               size = ggl.tg.sz, 
               colour = 'black',
               label = 'Google', 
               label.padding = unit(0.45, "lines"),
               hjust = 0, 
               vjust = 0) +
    geom_label(y = -Inf, x =  Inf, 
               size   = ggl.tg.sz, 
               colour = 'black',
               label.padding = unit(0.45, "lines"),
               label  = 'Imagery ©2022 TerraMetrics', 
               hjust  = 1, 
               vjust  = 0)  + 
    scaleBar(lon = box2[2] - scl.br.x.bf,      
             lat = box2[3] + scl.br.y.bf,      
             distanceLon = scl.dLon, 
             distanceLat = scl.dLat, 
             distanceLegend   = scl.dLeg, 
             dist.unit        = "km",
             legend.size      = scl.dLgsz,
             legend.colour    = "white",
             rec.colour       = "white",
             rec2.colour      = "white",
             arrow.length     = scl.arwLgt, 
             arrow.distance   = scl.arwDst, 
             arrow.North.size = scl.arwNsz,
             arrow.col = "white") +
    ggtitle(title)
  
 
  return(pp)
  
}