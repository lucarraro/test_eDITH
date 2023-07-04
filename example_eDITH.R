rm(list=ls())

library(rivnet)
library(terra)
library(LaplacesDemon)
library(BayesianTools)

river <- extract_river(outlet=c(637478,237413),
                    EPSG=21781, #CH1903/LV03 coordinate system
                    ext=c(6.2e5,6.6e5,2e5,2.5e5),
                    z=9)

river <- aggregate_river(river, maxReachLength = 2500)
hydrodata <- data.frame(data=c(8, 15), type=c("w","Q"), node=river$AG$outlet*c(1,1))
river <- hydro_river(hydrodata, river)

r1 <- rast(system.file("extdata/landcover.tif", package="rivnet"))
river <- covariate_river(r1, river)

ss <- sort(river$AG$A, index.return=T); ss <- ss$ix
samplingSites <- c(2,15,30,78,97,117,132,106,138,153,156,159,263,176,
  215,189,11,70,79,87,45,209,26,213) # don't touch!

tau <- 3*3600
log_p0 <- -16

covariates <- data.frame(urban=river$SC$locCov$landcover_1,
                         elev=river$AG$Z,
                         X=river$AG$X)
for (i in 1:length(covariates)) covariates[,i] <- (covariates[,i]-mean(covariates[,i]))/sd(covariates[,i])

param <- numeric(5)
param <- c(tau, log_p0, -2, -1, 0.5)
names(param) <- c("tau","log_p0","beta_urban","beta_Z","beta_X")

p <- eDITH:::eval.p(param, covariates)
C <- eDITH:::evalConc2_cpp(river, ss, river$AG$leng*river$AG$width, tau, p, "AG")

set.seed(1)
C_obs <- rnorm(length(samplingSites), C[samplingSites], 1e-13)

dd <- data.frame(ID=samplingSites, values=C_obs)

covariates <- data.frame(urban=river$SC$locCov$landcover_1,
                         agriculture=river$SC$locCov$landcover_2,
                         forest=river$SC$locCov$landcover_3,
                         elev=river$AG$Z,
                         log_drainageArea=log(river$AG$A))

if (!file.exists("out_BT.rda")){
out_BT <- eDITH::run_eDITH_BT(dd, river, covariates)
save(out_BT, file="out_BT.rda")
} else {load("out_BT.rda")}

out3 <- eDITH::run_eDITH_optim(dd, river, covariates)
save(out3, file="out_optim.rda")

out4 <- eDITH::run_eDITH_optim(dd, river)
save(out4, file="out_optim_noCov.rda")

out2 <- eDITH::run_eDITH_BT(dd, river) # works but doesn't converge

