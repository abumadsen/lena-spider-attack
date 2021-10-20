#30/06/20
#Mads F. Schou

#Lena comments:
# det gik lige op for mig at vi aldrig testede om der var forskel paa arterne i hvor lang tid det tog dem at angribe! Kan du koere saadan en test? Den burde vaere simplere fordi vi ikke behoever tids-aspekted som predictor. Efter et hurtigt plot i excel synes jeg ikke det ser ud som om der er meget forskel mellem arterne i deres respons til forskelligt str bytte, men det er en test der er vaerd at koere (jeg ved ikke lige hvorfor vi ikke har taenkt paa det noget foer!)
# 
# Det skal dog lige siges at der altsaa er en del angreb der aabenbart sker efter de standardiserede 10min, saa disse skal nok udelades fra denne analyse. 
# Jeg tror vi valgte at inkl. dem i hoved-analysen fordi vi jo der er mest interesserede i hvor mange gruppemedlemmer der angriber over tid, naar foerst byttet er bidt. (Men jeg skal vist lige inkludere lidt detaljer om det i metode-sektionen og ogsaa omregne mine overall acceptance rates). Nogle kolonier angriber indenfor 11min, og vi noterede dem naar vi kunne se aktivitet i reden omkring 9-10min og kunne se de nok snart ville angribe. Saa maaske prey acceptance skal oeges til 11min? Jeg kan faktisk ikke lige huske hvad vi gjorde med disse i DGAE paper, kan du? Saa vidt jeg umiddelbart kan se fra DGAE data saa kaldte vi et prey accepted op til 13 min, men jeg ved ikke om du aendrede paa det i analysen i R. 


pacman::p_load("dplyr","tidyr","doBy","MCMCglmm","parallel","coda","fitR")

DATAPATH = "interm"
OUTPATH = "code/5 Time to attack"

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

#Only one record per trial as we want to test time till first attack
mydat <- unique(mydat[, c("species","nestId","preysize","LatencyToAttackSec","trial","species2")])
mydat <- mydat[!is.na(mydat$LatencyToAttackSec),]
nrow(mydat) #192

#Confirm that only one obsevation per trial
mydat$trial.nest <- paste(mydat$trial,mydat$nestId, sep = "_")
unique(table(mydat$trial.nest)) #1 -> confirmed!
table(table(mydat$nestId)) #But some nests are tested multiple times, so need this as a random effect

#Removing time points after 10 min da det ikke længere kan ses som et attack
mydat <- mydat[mydat$LatencyToAttackSec < 10*60,]
nrow(mydat) #178



####################################
##---- Prepping (Keep minor at this stage)
####################################

mydat <- mydat %>%
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate_at(c("species","species2","nestId","trial"), list(~factor(.)))


####################################
##---- TEST RUN
####################################
# # 
# MyStart <- Sys.time()
# mtest <- MCMCglmm(LatencyToAttackSec ~ species-1 + species:preysize_z,
#                   random = ~ nestId,
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


m1_3.2way.nestid <- mclapply(1:3, function(i) {
  MCMCglmm(LatencyToAttackSec ~ species-1 + species:preysize_z,
           random = ~  nestId,
           data   = mydat,
           family = "poisson",
           prior  = MyPrior,
           #thin   = 1,burnin = 0,nitt   = 1000,
           nitt=2500000, thin=2500, burnin=100000,
           verbose = T)
}, mc.cores=3)

m1_3.1way.nestid <- mclapply(1:3, function(i) {
  MCMCglmm(LatencyToAttackSec ~ species-1 + preysize_z,
           random = ~ nestId,
           data   = mydat,
           family = "poisson",
           prior  = MyPrior,
           #thin   = 1,burnin = 0,nitt   = 1000,
           nitt=2500000, thin=2500, burnin=100000,
           verbose = T)
}, mc.cores=3)


save.image(paste(OUTPATH,"/", "Time to attack.RData", sep = ""))
