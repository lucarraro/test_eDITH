rm(list=ls())
library(rivnet)
library(OCNet)
library(Rcpp)
library(terra)
library(BayesianTools)
library(LaplacesDemon)
library(readxl)
library(fields)
library(eDITH) # from github.com/lucarraro/eDITH

#sourceCpp("evalConc2.cpp")
#source("run_eDITH_func.R")

if (!dir.exists("results")) dir.create("results")

if (!file.exists("data/Thur.rda")){
  Thur <- rivnet::extract_river(outlet = c(735010, 261530),
                                EPSG = 21781,
                                ext = c(700000, 770000, 220000, 270000),
                                z = 11)

  Thur <- rivnet::aggregate_river(Thur, thrA = 0.25e6, maxReachLength = 1000)

  hydroData <- data.frame(x=c(723675, 727110, 718840, 737270),
                          y=c(252720, 247290, 248440, 251290),
                          w=c(    35,     20,    2.5,      8),
                          d=c(1.2034, 0.6784,     NA, 0.2500),
                          Q=c(52.8667,7.9083, 0.2595, 1.3392))

  siteAG <- numeric(length(hydroData$x))
  for (i in 1:length(hydroData$x)){
    ss <- locate_site(hydroData$x[i], hydroData$y[i], Thur)
    siteAG[i] <- ss$AGnode
    #readline(prompt="Press [enter] to continue")
  }
  hyData <- data.frame(type=c(rep("w",4),rep("d",3),rep("Q",4)),
                       data=c(hydroData$w, hydroData$d[!is.na(hydroData$d)], hydroData$Q),
                       node=c(siteAG, siteAG[!is.na(hydroData$d)], siteAG))
  Thur <- hydro_river(hyData, Thur)

  #r1 <- terra::rast("data/landcover.asc")

  landCovCH <- read.csv("C:/Users/carrarlu/Documents/ag-b-00.03-37-area-csv.csv", sep = ";") # add this .csv file to the source folder
  landCovCH <- terra::rast(data.frame(landCovCH$E, landCovCH$N, landCovCH$LU18_4),
                           type = "xyz", crs = "EPSG:2056") # convert into raster
  # Legend: 1-Urban; 2-Agriculture; 3-Forest; 4-Improductive
  landCovCH <- terra::project(landCovCH, crs("EPSG:21781"))
  Thur <- rivnet::covariate_river(landCovCH, Thur, overwrite=TRUE)
  rm(landCovCH)

  geol_rast <- terra::rast("data/Clipped_Geology.asc")
  crs(geol_rast) <- "EPSG:21781"
  vv <- values(geol_rast)
  vv[vv==7 | vv==5 | vv==2] <- 1001 #"alluvial"
  vv[vv==30 | vv==57 | vv==108 | vv==48 | vv==45] <- 1002 # "alpine"
  vv[vv==14] <- 1003 # "loess"
  vv[vv==1 | vv==3 | vv==8 | vv==9 |  vv==20 |  vv==27] <- 1004 #"molasses"
  vv[vv==6] <- 1005 #"moraines"
  vv[vv==37 | vv==40 | vv==62 | vv==69 | vv==76 | vv==87] <- 1006 # "other"
  vv[vv==4] <- 1007 #"peat"
  vv[vv==35 | vv==43 | vv==97] <- 1008 #"scree"
  vv[vv==106] <- 1009 # "water"
  values(geol_rast) <- vv

  Thur <- rivnet::covariate_river(geol_rast, Thur)

  names(Thur$SC$locCov) <- names(Thur$SC$upsCov) <-c("urban", "agriculture", "forest", "improductive",
                                                     "alluvial","alpine","loess","molasses","moraines",
                                                     "other","peat","scree","water")

  save(Thur,file="Thur.rda",compress="xz")
} else {load("data/Thur.rda")}

# geographical covariates
siteClusters <- read_excel("data/siteClusters.xlsx")
geo_cluster <- numeric(Thur$AG$nNodes)
geoCov <- data.frame(matrix(0,Thur$AG$nNodes,dim(siteClusters)[1]))
for (i in 1:dim(siteClusters)[1]){
  geo_cluster[Thur$AG$upstream[[siteClusters$AG[i]]]] <- i
}
for (i in 1:dim(siteClusters)[1]){
  geoCov[geo_cluster==i, i] <- 1
}
names(geoCov) <- siteClusters$ID

covariates <- cbind(data.frame(logDrainageArea=log(Thur$AG$A), # streamOrder=Thur$AG$streamOrder,
                               elevation=Thur$AG$Z, slope=Thur$AG$slope),
                    Thur$SC$locCov[c("urban","agriculture","forest")],
                    Thur$SC$upsCov[c("alluvial","alpine","molasses","moraines","peat")], #"water", "loess", "scree"
                    geoCov)

# sampling sites
samplingSites <- read_excel("data/Coordinates_new.xlsx")
X <- samplingSites$X_new; Y <- samplingSites$Y_new
sampSiteAG <- length(X)
X[4] <- 726000;
X[7] <- 723700; Y[7] <- 254100
X[10] <- 733350
Y[15] <- 256300
X[21] <- 728400; Y[21] <- 251900
X[27] <- 737700; Y[27] <- 247200
X[39] <- 731300
Y[43] <- 239100
for (i in 1:length(X)){
  tmp <- locate_site(X[i], Y[i], Thur,showPlot=F)
  #title(sprintf("Site %d  -  X %d  -  Y %d",samplingSites$SiteID[i], X[i], Y[i]))
  sampSiteAG[i] <- tmp$AGnode
  #readline("Press Enter to continue:")
}
X_site <- X; Y_site <- Y;
site_key <- sort(samplingSites$SiteName, index.return=T); site_key <- site_key$ix

allData_nP <- read.csv("data/MZBTaxa_summed_notPooled.csv") # not pooled
familyNames <- as.vector(allData_nP["Family"]); familyNames <- familyNames$Family
allData_nP <- allData_nP[ ,-1]
row.names(allData_nP) <- familyNames
allData_nP <- as.data.frame(t(allData_nP))
allData_nP <- allData_nP[-which(row.names(allData_nP)=="V26.5"), ] # remove additional run on V26
rn <- row.names(allData_nP); rn[rn=="H27.3"] <- "L27.3" # correct mislabeled run
row.names(allData_nP) <- rn


kols <- hcl.colors(1000,rev=T)


for (ff in familyNames){
  cat(sprintf('%d) %s  -  %s \n',which(familyNames==ff),ff,date()))
  if(!file.exists(paste0("results/",ff,".rda"))){
    out <- NULL
    eval(parse(text=paste0('save(out, file="results/',ff,'.rda")')))

    data <- data.frame(ID=sampSiteAG[match(substr(row.names(allData_nP),1,3), samplingSites$RB_Site_code)],
                       values=allData_nP[[ff]]) # NON-POOLED DATA

    out <- run_eDITH_BT(data, Thur, covariates, no.det=FALSE, ll.type = "nbinom",
                        mcmc.settings = list(iterations = 2.7e4, burnin=1.8e4, message = TRUE, thin = 10)) # short version
    out[["taxon"]] <- ff
    out <- eDITH::eval_posterior_eDITH(out, Thur, covariates)

    pdf(paste0("results/",ff,"_map.pdf"),width=20/2.54,height=20/2.54)
    plot(out$probDetection_median, Thur, colLevels=c(0,1)); title(ff)
    for(i in 1:length(unique(data$ID))){
      value <- sum(data$values[data$ID==unique(data$ID)[i]])
      if (value>0){
        ind_kol <- floor(log10(value)/4*999)+1
        if (ind_kol>1000) ind_kol <- 1000
        points(Thur$AG$X[unique(data$ID)[i]], Thur$AG$Y[unique(data$ID)[i]], pch=21, bg=kols[ind_kol])
      } else {points(Thur$AG$X[unique(data$ID)[i]], Thur$AG$Y[unique(data$ID)[i]], pch=21,  bg="white")}
    }
    imagePlot(col=kols, zlim=c(0,4),  legend.only=T,smallplot=c(0.2,0.8,0.1,0.12),
              horizontal=T, legend.lab="No. pooled reads (log10)")
    dev.off()

    pdf(paste0("results/",ff,"_LP.pdf"),width=18/2.54,height=12/2.54)
    plot(out$outMCMC$chain[[1]][-1,"LP"],type="l",ylab="Log-posterior"); title(ff)
    lines(out$outMCMC$chain[[2]][-1,"LP"],col="red",type="l")
    lines(out$outMCMC$chain[[3]][-1,"LP"],col="blue",type="l")
    dev.off()

    eval(parse(text=paste0('save(out, file="results/',ff,'.rda", compress="xz")')))

    rm(ff)
    gc(verbose=F)
  }
}


# for (i in 1:dim(covariates)[2]){
#   nam <- names(covariates)[i]
#   fnam <- paste0("covariate_plot/",nam,".pdf")
#    pdf(file=fnam,width=20/2.54,height=20/2.54)
#    plot(Thur,covariates[,i]); title(nam)
#    dev.off()
# }


