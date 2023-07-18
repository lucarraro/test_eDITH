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

#try with more sites
#samplingSites <- 1:river$AG$nNodes

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
C_obs <- rnorm(3*length(samplingSites), rep(C[samplingSites],3), 3e-13)
C_obs[C_obs<8e-13] <- 0

Read_obs <- rgeom(3*length(samplingSites), prob=rep(1/(1+1e11*C[samplingSites]), 3))
omega <- 20
Read_obs <- rnbinom(3*length(samplingSites),
                    size = rep(1e13*C[samplingSites]/(omega-1),3), prob = 1/omega)

dd <- data.frame(ID=rep(samplingSites,3), values=C_obs)
dd.read <- data.frame(ID=rep(samplingSites,3), values=Read_obs)

covariates <- data.frame(urban=river$SC$locCov$landcover_1,
                         agriculture=river$SC$locCov$landcover_2,
                         forest=river$SC$locCov$landcover_3,
                         elev=river$AG$Z,
                         log_drainageArea=log(river$AG$A))

if (!file.exists("out_BT.rda")){
out_BT <- eDITH::run_eDITH_BT(dd, river, covariates)
save(out_BT, file="out_BT.rda")
} else {load("out_BT.rda")}

out2 <- eDITH::run_eDITH_BT(dd, river, covariates, no.det=TRUE)

out.read <- eDITH::run_eDITH_BT(dd.read, river, covariates, ll.type="nbinom")
out.geom <- eDITH::run_eDITH_BT(dd.read, river, covariates, ll.type="geom")
out.geom.noDet <- eDITH::run_eDITH_BT(dd.read, river, covariates, ll.type="geom", no.det=T)

out.geom.short <- eDITH::run_eDITH_BT(dd.read, river, covariates, ll.type="geom",
                                      mcmc.settings = list(iterations = 3e3, message = T, thin = 10))

## alternative functions and options
out3 <- eDITH::run_eDITH_optim(dd, river, covariates) # this function finds a single best-fit parameter set

# use AEMs as covariates
out4 <- eDITH::run_eDITH_optim(dd, river)
out5 <- eDITH::run_eDITH_BT(dd, river)

# append the first 6 AEMs to the provided covariates
out6 <- eDITH::run_eDITH_BT(dd, river, covariates, use.AEM=TRUE, n.AEM=6)

# produce posterior sample to be used as new prior
set.seed(1)
outSample <- eDITH::run_eDITH_BT(dd, river, covariates,
                              mcmc.settings=list(iterations=9e5, burnin = 6e5, message = TRUE, thin = 30))
save(outSample, file="outSample.rda", compress="xz")

pp <- createPriorDensity(outSample$outMCMC)
names(pp$lower) <- names(pp$upper) <- colnames(outSample$outMCMC$chain[[1]])[1:8]
# the three last columns are for log-posterior, log-likelihood, log-prior
out.new <- run_eDITH_BT(dd, river, covariates, prior=pp)
