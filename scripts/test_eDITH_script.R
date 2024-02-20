#-----------------------------------------------------------------------------#

# Illustrating the eDITH package with two case studies

# Script accompanying the manuscript "eDITH: leveraging environmental DNA to
#   spatially project biodiversity and ecological indices across river networks"
#   Carraro & Altermatt, 2024

# author: Luca Carraro (luca.carraro@eawag.ch)
# University of Zurich & Eawag, Switzerland
# date: January 5, 2024

#-----------------------------------------------------------------------------#

rm(list=ls())

library(rivnet)
library(eDITH)
library(DHARMa)
library(BayesianTools)
library(terra)

riverData <- read.csv(file="../data/riverData.csv") # file containing location and hydrological data for the two rivers
thr_species <- 5 # exclude species with less then thr_species detections from modelling

# initialize variables
probDet_median <- probDet_025 <- probDet_975 <- vector("list",2)
signifCovariates <- gelmanDiag <- SST <- speciesID <- vector("list",2)

names(probDet_median) <- names(probDet_025) <- names(probDet_975) <- riverData$river
names(signifCovariates) <- names(gelmanDiag) <- names(SST) <- names(speciesID) <- riverData$river

n_covariates <- 10 # number of covariates used in eDITH

for (i in 1:length(riverData$river)){ # loop on rivers
  ID <- riverData$river[i]
  message(paste0('river ',ID,'\n'),appendLF=F)
  if (!file.exists(paste0('../data/river',ID,'.rda'))){ # if river object was already generated, skip
    # otherwise, extract river
    river <- rivnet::extract_river(outlet = c(riverData$X.outlet[i], riverData$Y.outlet[i]),
                                   EPSG = riverData$EPSG[i],
                                   ext = c(riverData$X.min[i], riverData$X.max[i],
                                           riverData$Y.min[i], riverData$Y.max[i]),
                                   z = riverData$z[i],
                                   showPlot = T,
                                   threshold_parameter = 1000,
                                   displayUpdates = 1,
                                   n_processes = 8)
    # aggregate river
    river <- rivnet::aggregate_river(river, thrA=riverData$thrA[i],
                                     maxReachLength=1000,
                                     equalizeLengths=TRUE)
    # create data frame with hydrological variables
    hd <- data.frame(data=c(riverData$Q.outlet[i], riverData$w.outlet[i]),
                     type=c("Q", "w"),
                     node=c(1,1)*river$AG$outlet)
    # extrapolate hydrological variables
    river <- hydro_river(hd, river)

    # identify AG nodes where sampling sites are located
    samplingSites <- read.csv(file=paste0("../data/samplingSites",ID,".csv"))
    if (i==1) samplingSites$X[10] <- samplingSites$X[10] - 200 # fix position of site 10 for catchment B
                                                               # otherwise it would be assigned to the wrong reach
    AGnode <- numeric(length(samplingSites$X))
    for (j in 1: length(samplingSites$X)){
      tmp <- locate_site(samplingSites$X[j], samplingSites$Y[j], river)#, showPlot = T)
      #title(samplingSites$siteID[j])
      #readline('Press Enter to continue: ')
      AGnode[j] <- tmp$AGnode
    }
    samplingSites[["AGnode"]] <- AGnode

    # save river object and sampling sites
    eval(parse(text=paste0('river',ID,' <- river')))
    eval(parse(text=paste0('samplingSites',ID,' <- samplingSites')))
    eval(parse(text=paste0('save(river',ID,',samplingSites',ID,',file="../data/river',ID,'.rda")' )))

  } else { # load existing river object if already generated
    eval(parse(text=paste0('load("../data/river',ID,'.rda")' )))
    eval(parse(text=paste0('river <- river',ID)))
    eval(parse(text=paste0('samplingSites <- samplingSites',ID)))
  }

  # read species-by-site table
  SST[[i]] <- read.csv(paste0("../data/siteSpeciesTable",ID,".csv"))
  SST[[i]][is.na(SST[[i]])] <- 0
  speciesID[[i]] <- which(rowSums(SST[[i]][,-1]!=0) >= thr_species) # choose species with at least thr_species presences

  # initialize variables for export
  probDet_median[[i]] <- data.frame(matrix(0,river$AG$nNodes, length(speciesID[[i]])))
  probDet_025[[i]] <- probDet_975[[i]] <- probDet_median[[i]]
  names(probDet_025[[i]]) <- names(probDet_975[[i]]) <- names(probDet_median[[i]]) <-SST[[i]][speciesID[[i]],1]

  signifCovariates[[i]] <- data.frame(matrix(0,n_covariates,length(speciesID[[i]])))
  names(signifCovariates[[i]]) <- SST[[i]][speciesID[[i]],1]

  row.names(signifCovariates[[i]]) <- paste0("AEM",1:n_covariates)
  gelmanDiag[[i]] <- data.frame(matrix(0,n_covariates+3,length(speciesID[[i]])))

  names(gelmanDiag[[i]]) <- SST[[i]][speciesID[[i]],1]
  row.names(gelmanDiag[[i]]) <- c("tau","log_p0",paste0("AEM",1:n_covariates),"omega")

  for (iS in 1:length(speciesID[[i]])){ # loop over species
    indSp <- speciesID[[i]][iS]
    nam <- SST[[i]][indSp,1]
    fnam <- paste0('../results',ID,'/',nam,'.rda')
    message(paste0(' ',nam,'\n'),appendLF=F)
    if(!file.exists(fnam)){ # if eDITH results were already generated, skip
      out <- NULL
      save(out,file=fnam) # save dummy file (enables parallel computing)
      values <- as.numeric(SST[[i]][indSp,-1])
      data <- data.frame(ID=samplingSites$AGnode, values=values) # create data frame with eDNA data
                                                                 # for the species at hand
      # run eDITH model via BayesianTools
      out <- run_eDITH_BT(data, river, ll.type="nbinom",
                          n.AEM=n_covariates, par.AEM=list(weight="exponential"),
                          verbose=T)
      # evaluate quantiles from posterior distribution
      out <- eval_posterior_eDITH(out, river, quant=c(0.025,0.5, 0.975))
      # run posterior predictive simulations
      pps <- posterior_pred_sim_eDITH(out, river)
      out$outMCMC <- NULL # delete runMCMC output to save disk space
      save(out,pps,file=fnam) # save file
    } else {load(fnam)    } # load eDITH output if already generated
    message('\n',appendLF=F)

    fnam_pdf <- paste0('../results',ID,'/',nam,'.pdf')
    fnam_dharma <- paste0('../results',ID,'/',nam,'_dharma.pdf')

    if(!is.null(out$probDetection_quantile)){
      # store posterior detection probability (median, 0.025- and 0.975-quantiles),
      #       Gelman diagnostics and significance of covariates
      probDet_median[[i]][,iS] <- out$probDetection_quantile[2,]
      probDet_025[[i]][,iS] <- out$probDetection_quantile[1,]
      probDet_975[[i]][,iS] <- out$probDetection_quantile[3,]
      gelmanDiag[[i]][,iS] <- out$gD$psrf[,1]
      for (j in 1:length(out$covariates)){
        if (out$cI[1,j+2]>0) signifCovariates[[i]][j,iS] <- 1
        if (out$cI[2,j+2]<0) signifCovariates[[i]][j,iS] <- -1
      }

      # plot maps of modelled detection probability and read counts
      if (!file.exists(fnam_pdf)){
        pdf(file=fnam_pdf, width=30/2.54, height=20/2.54)
        par(mfrow=c(1,2))
        plot(out$probDetection_quantile[2,],  river,
             args_imagePlot=list(legend.lab="Posterior median detection probability"))
        points_colorscale(samplingSites$X, samplingSites$Y, log10(out$data$values),
                          cex=0.75, force.range = FALSE,
                          bg.palette=hcl.colors(1000, "Blues 3", rev=T),
                          horizontal=T, legend.lab = "Observed read counts (log 10)")
        title(nam)

        plot(out$C_quantile[2,], river, colPalette=hcl.colors(1000, "Reds 3",rev=T),
             addLegend=F)
        points_colorscale(samplingSites$X, samplingSites$Y, out$data$values, cex=0.75,
                          bg.range=range(out$C_quantile[2,]), horizontal=T,
                          legend.lab = "Read counts (observed vs. posterior median)")
        title(nam)
        dev.off()
      }
      # plot DHARMa diagnostics
      if (!file.exists(fnam_dharma)){
        pdf(file=fnam_dharma, width=30/2.54, height=20/2.54)
        out.sim <- createDHARMa(pps, out$data$values)
        plot(out.sim)
        dev.off()
      }
    }
  }
}

# plot alpha diversity (Fig. 4 of the manuscript)
pdf("Fig4.pdf",width=30/2.54,height=25/2.54)
sp_name <- c("Cottus bairdii","Eleotris oxycephala") #c(28, 27) #
deu <- colorRampPalette(c("yellow","red","black"))
north_length <- c(5000,1000)
par(mfrow=c(2,2))
for (i in 1:length(riverData$river)) {
  ID <- riverData$river[i]
  eval(parse(text=paste0('river <- river',ID)))
  eval(parse(text=paste0('samplingSites <- samplingSites',ID)))
  plot(probDet_median[[i]][[sp_name[i]]],  river,
       colPalette=hcl.colors(1000, "Reds 3", rev=T), colLevels = c(0,1),
       args_imagePlot=list(legend.lab="Posterior median detection probability"))
  values <- SST[[i]][which(SST[[i]][,1]==sp_name[i]),-1]
  points_colorscale(samplingSites$X, samplingSites$Y, log10(as.numeric(values)),
                    cex=1, force.range = FALSE,
                    bg.palette=hcl.colors(1000, "Blues 3", rev=T),
                    bg.range=c(1,5),
                    horizontal=T, legend.lab = "Observed read counts (log 10)")
  title(sp_name[i])

  plot(rowSums(probDet_median[[i]]>0.5), river, addLegend=F, max_lwd=3, colLevels=c(0,25))
  terra::sbar(d = 1000,
              xy = c(min(river$FD$X)+0.05*(max(river$FD$X)-min(river$FD$X)),
                     min(river$FD$Y)+0.7*(max(river$FD$Y)-min(river$FD$Y))),
              label = "1 km")
  terra::north(d = north_length[i],
               xy = c(min(river$FD$X)+0.05*(max(river$FD$X)-min(river$FD$X)),
                      min(river$FD$Y)+0.9*(max(river$FD$Y)-min(river$FD$Y))))
  rivnet::points_colorscale(samplingSites$X, samplingSites$Y, colSums(SST[[i]][speciesID[[i]],-1]>0),
                    bg.palette = deu(1000),bg.range=c(0,25),
                    cex=1, horizontal = T, legend.lab = "Species richness")
}
dev.off()

# Plot range of 95%-equal tailed credible interval for detection probability
par(mfrow=c(1,2))
plot(probDet_025[[1]][[sp_name[1]]],  riverSydenham,
      colPalette=hcl.colors(1000, "Reds 3", rev=T), colLevels = c(0,1),
      addLegend = FALSE)
plot(probDet_975[[1]][[sp_name[1]]],  riverSydenham,
     colPalette=hcl.colors(1000, "Reds 3", rev=T), colLevels = c(0,1),
     args_imagePlot=list(legend.lab="C. bairdii - Posterior detection probability"))

