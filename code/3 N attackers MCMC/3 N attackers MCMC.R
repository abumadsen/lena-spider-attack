#30/06/20
#Mads F. Schou

rm(list = ls(all = TRUE))
pacman::p_load("dplyr","tidyr","pedantics","doBy","MCMCglmm","parallel","coda","fitR")

DATAPATH = "Temp"
OUTPATH = "Analyses/3 N attackers MCMC"

#---- Data
mydat <- read.table(paste(DATAPATH,"social attack data for Mads 22.05.20_prepped.csv",sep = "/"), sep = ",", header = TRUE)

#---- Prior
MyPrior <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),
  G=list(G1=list(V        = diag(2),
                 nu        = 1.002)))

####################################
##---- Filtering (Keep minor at this stage)
####################################

#Two highly influenctial
outlier_ids <- mydat$trial[mydat$attackers > 20]
mydat <- mydat[!mydat$trial %in% outlier_ids,]

####################################
##---- Prepping (Keep minor at this stage)
####################################

mydat <- mydat %>%
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate(timepoint_z = scale(timepoint)[,1]) %>%
  mutate_at(c("species","species2","nestId","trial"), list(~factor(.))) %>%
  mutate(preysize.log = log(preysize)) %>%
  mutate(preysize.log_z = scale(preysize.log)[,1])
  
  

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

m1.3way <- MCMCglmm(attackers ~ species-1 + species:timepoint_z+ species:preysize_z+species:preysize_z:timepoint_z,
               random = ~ us(1+timepoint_z):trial,
               data   = mydat,
               family = "poisson",
               prior  = MyPrior,
               #thin   = 1,burnin = 0,nitt   = 1000,
               nitt=2500000, thin=2500, burnin=100000,
               verbose = T)


m1.3way.log <- MCMCglmm(attackers ~ species-1 + species:timepoint_z+ species:preysize.log_z+species:preysize.log_z:timepoint_z,
                    random = ~ us(1+timepoint_z):trial,
                    data   = mydat,
                    family = "poisson",
                    prior  = MyPrior,
                    #thin   = 1,burnin = 0,nitt   = 1000,
                    nitt=2500000, thin=2500, burnin=100000,
                    verbose = T)

m1.2way <- MCMCglmm(attackers ~species-1 + species:timepoint_z+ species:preysize.log_z+preysize.log_z:timepoint_z,
                    random = ~ us(1+timepoint_z):trial,
                    data   = mydat,
                    family = "poisson",
                    prior  = MyPrior,
                    #thin   = 1,burnin = 0,nitt   = 1000,
                    nitt=2500000, thin=2500, burnin=100000,
                    verbose = T)
# 
# m1.1way <- MCMCglmm(attackers ~ species-1 + preysize_z + timepoint_z,
#                     random = ~ us(1+timepoint_z):trial,
#                     data   = mydat,
#                     family = "poisson",
#                     prior  = MyPrior,
#                     #thin   = 1,burnin = 0,nitt   = 1000,
#                     nitt=530000, thin=500, burnin=30000,
#                     verbose = T)



####################################
##---- FINAL RUN FOR PLOTTING
####################################
#This is necessary as I cannot get predict to run with the random slopes

MyPriorPlot <- list(
  R=list(V = diag(1), nu=0.002, fix = FALSE),
  G=list(G1=list(V        = diag(1),
                 nu        = 0.002)))
m1.3wayPlot <- MCMCglmm(attackers ~ species-1 + species:timepoint_z+ species:preysize_z+species:preysize_z:timepoint_z,
                    random = ~ trial,
                    data   = mydat,
                    family = "poisson",
                    prior  = MyPriorPlot,
                    #thin   = 1,burnin = 0,nitt   = 1000,
                    nitt=2500000, thin=2500, burnin=100000,
                    verbose = T)

m1.3way.logPlot <- MCMCglmm(attackers ~ species-1 + species:timepoint_z+ species:preysize.log_z+species:preysize.log_z:timepoint_z,
                        random = ~ trial,
                        data   = mydat,
                        family = "poisson",
                        prior  = MyPriorPlot,
                        #thin   = 1,burnin = 0,nitt   = 1000,
                        nitt=2500000, thin=2500, burnin=100000,
                        verbose = T)

save.image(paste(OUTPATH,"/", "N attackers no outliers.RData", sep = ""))
