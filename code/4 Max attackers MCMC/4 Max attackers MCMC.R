#30/06/20
#Mads F. Schou

rm(list = ls(all = TRUE))
pacman::p_load("dplyr","tidyr","pedantics","doBy","MCMCglmm","parallel","coda","fitR")

DATAPATH = "Temp"
OUTPATH = "Analyses/4 Max attackers MCMC"

#---- Data
mydat <- read.table(paste(DATAPATH,"social attack data for Mads 22.05.20_prepped.csv",sep = "/"), sep = ",", header = TRUE)

#---- Prior
MyPrior <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),
  G=list(G1=list(V        = diag(1),
                 nu        = 0.002)))

####################################
##---- Filtering (Keep minor at this stage)
####################################

#Two highly influenctial
outlier_ids <- mydat$trial[mydat$attackers > 20]
mydat <- mydat[!mydat$trial %in% outlier_ids,]

####################################
##---- Prepping (Keep minor at this stage)
####################################

# Note that we could leave nestID out here, as most only have one observation per nest id.
mydat <- summaryBy(attackers ~ species + nestId + preysize + trial + preysize_z + species2, data = mydat, FUN = max)

mydat <- mydat %>%
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate_at(c("species","nestId","trial"), funs(factor(.)))

####################################
##---- TEST RUN
####################################
# 
# MyStart <- Sys.time()
# mtest <- MCMCglmm(attackers ~ timepoint_z + species-1 + species:preysize_z,
#                   random = ~ us(1+timepoint_z):trial,
#                   data   = mydat,
#                   family = "poisson",
#                   prior  = MyPrior,
#                   nitt=100000, thin=1, burnin=0,
#                   verbose = T)
# TotalTime <- Sys.time() - MyStart
# 
# #---------- Inspect how we need to set thin and burnin to get efficient sampling
# 
# #Extract trace
# Mytrace <- mtest$VCV
# mcmc.trace <- mcmc(Mytrace)
# 
# #Plot trace plots
# plot(mcmc.trace)
# 
# #Cut by burnin
# mcmc.trace.burned <- burnAndThin(mcmc.trace, burn = 30000)
# plot(mcmc.trace.burned) #plot again
# 
# #Autocorr
# autocorr.plot(mcmc.trace.burned)
# autoOut <- autocorr.diag(mcmc.trace.burned)
# autoOut <- data.frame("rowMax" = apply(autoOut,1,max),autoOut) #Get max auto across all variates
# autoOut
# 
# #Thin it
# mcmc.trace.burned.thinned <- burnAndThin(mcmc.trace.burned, thin = 500)
# autocorr.plot(mcmc.trace.burned.thinned) #Inspect autocorr again
# autoOut <- autocorr.diag(mcmc.trace.burned.thinned)
# autoOut <- data.frame("rowMax" = apply(autoOut,1,max),autoOut) #Get max auto across all variates
# autoOut
# 
# plot(mcmc.trace.burned.thinned) # Looks good

####################################
##---- FINAL RUN
####################################

#---------------------------- No poly

MyStart <- Sys.time()
#set.seed(1)
#m1_3 <- mclapply(1:3, function(i) {
m1 <- MCMCglmm(attackers.max ~ species-1 + species:preysize_z,
               random = ~ nestId,
               data   = mydat,
               family = "poisson",
               prior  = MyPrior,
               #thin   = 1,burnin = 0,nitt   = 1000,
               nitt=1030000, thin=500, burnin=30000,
               verbose = T)
#}, mc.cores=3)
TotalTime <- Sys.time() - MyStart 


#---------------------------- Poly

#set.seed(1)
#m1_3 <- mclapply(1:3, function(i) {
m1.poly <- MCMCglmm(attackers.max ~ species-1 + species:poly(preysize_z,2),
               random = ~ nestId,
               data   = mydat,
               family = "poisson",
               prior  = MyPrior,
               #thin   = 1,burnin = 0,nitt   = 1000,
               nitt=1030000, thin=500, burnin=30000,
               verbose = T)
#}, mc.cores=3)



save.image(paste(OUTPATH,"/", "Max attackers no outliers.RData", sep = ""))

