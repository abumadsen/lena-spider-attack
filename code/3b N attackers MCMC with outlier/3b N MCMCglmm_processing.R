#30/06/20
#Mads F. Schou
#Example of data processing of output from MCMCglmm

pacman::p_load("reshape","openxlsx","MCMCglmm")
source("https://raw.githubusercontent.com/abumadsen/custom-R-functions/main/MCMCglmm_processing_functions.R")

#1. Model Diagnostics
#2. Report Results

PathToModel = "code/3 N attackers MCMC/"
FileName = "N attackers no outliers.nest.RData"
#Load run
load(paste(PathToModel,FileName, sep = ""))
PathOutput = "results/3 N attackers MCMC/with-nest/"

#If one:
#m1 <- m1.3way.nest

#If multiple chains, isolate one:
m1 <- m1_3.3way.nestid[[1]]
m1_3 <- m1_3.3way.nestid

#obtain stats
burnin <- summary(m1)$cstats[1]-1
thin <- summary(m1)$cstats[3]
itt <- summary(m1)$cstats[2]+thin-1
specs <- paste("",thin,burnin,itt, sep = "_")

#######################################
###--- 1. Model Diagnostics
#######################################

#---------- Check autocorrelation
autoOut <- autocorr.diag(m1$VCV)
autoOut <- data.frame("rowMax" = apply(autoOut,1,max),autoOut) #Get max auto across all variates
write.csv(autoOut, paste(PathOutput, "autocorrelation",specs,".csv", sep = ""), quote = FALSE, row.names = TRUE)
#autocorr.plot(mcmc(m1$VCV),lag.max = 5)

#-------------- if multiple chains

# --- Gelman Rubin Criterion
# Fixed effects
m1_3Sol <- lapply(m1_3, function(m) m$Sol[,colnames(m$Sol) %in% row.names(summary(m)$solutions)]) #if pr=TRUE this removes all unwanted posteriors)
m1_3Sol <- do.call(mcmc.list, m1_3Sol)
gelman.diag(m1_3Sol)
pdf(paste(PathOutput,"GRC_fixed",specs,".pdf",sep = ""), family = "Times")
#par(mfrow=c(4,2), mar=c(2,2,1,2));gelman.plot(m1_3Sol, auto.layout=F,autoburnin = FALSE)
par(mfrow=c(4,2), mar=c(2, 1, 1, 1)); plot(m1_3Sol, ask=F, auto.layout=F)
dev.off()

# Random effects
m1_3VCV <- lapply(m1_3, function(m) m$VCV)
m1_3VCV <- do.call(mcmc.list, m1_3VCV)
gelman.diag(m1_3VCV, multivariate = FALSE)
pdf(paste(PathOutput,"GRC_random",specs,".pdf",sep = ""), family = "Times")
#par(mfrow=c(4,2), mar=c(2,2,1,2))
#gelman.plot(m1_3VCV, auto.layout=F)
par(mfrow=c(4,2), mar=c(2, 1, 1, 1))
plot(m1_3VCV, ask=F, auto.layout=F)
dev.off()

#-------------- if one chain

# #--- Plot chain of random effects
# pdf(paste(PathOutput,"chains_random",specs,".pdf", sep = ""), family = "Times")
# par(mfrow=c(5,2), mar=c(2,2,1,0))
# plot(m1$VCV, auto.layout=F)
# dev.off()
# 
# #--- Plot chain of fixed effects
# pdf(paste(PathOutput,"chains_fixed",specs,".pdf", sep = ""), family = "Times")
# par(mfrow=c(5,2), mar=c(2,2,1,0))
# plot(m1$Sol[,colnames(m1$Sol) %in% row.names(summary(m1)$solutions)], auto.layout=F)
# dev.off()


#######################################
###--- 2. Report Results
#######################################

#Open workbook
OpenWorkBookMCMCglmm(title = "Table S1: Results from analysis of the number of attackers", sheetname = "Results")
#Add results from Fixed, Random and Correlations
MyModel <- m1
FixedOut <- ReportFixedMCMC(MyModel, remove = "")
RandomOut <- ReportRandomVarianceMCMC(MyModel)
CorrOut <- ReportCorrelationsMCMC(MyModel)
writeData(workbook, "Results", FixedOut, startCol = 1, startRow = 2,headerStyle = hs2,colNames =TRUE)
writeData(workbook, "Results", RandomOut, startCol = 1, startRow = nrow(FixedOut)+3,headerStyle = hs2)
writeData(workbook, "Results", CorrOut, startCol = 1, startRow = nrow(FixedOut)+nrow(RandomOut)+4,headerStyle = hs2)

#Testing for differences in change in number of attackers with increasing preysize among the three species
MyCont1 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesdumicola:preysize_z", Eff2 = "speciessarasinorum:preysize_z")
MyCont2 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesmimosarum:preysize_z", Eff2 = "speciessarasinorum:preysize_z")
MyCont3 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesmimosarum:preysize_z", Eff2 = "speciesdumicola:preysize_z")
MyCont <- rbind(MyCont1,MyCont2,MyCont3)
writeData(workbook, "Results", MyCont, startCol = 1, startRow = nrow(FixedOut)+nrow(RandomOut)+nrow(CorrOut)+5,headerStyle = hs2)

#Testing for differences in change in number of attackers with increasing time among the three species
MyCont1 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesdumicola:timepoint_z", Eff2 = "speciessarasinorum:timepoint_z")
MyCont2 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesmimosarum:timepoint_z", Eff2 = "speciessarasinorum:timepoint_z")
MyCont3 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesmimosarum:timepoint_z", Eff2 = "speciesdumicola:timepoint_z")
MyContTime <- rbind(MyCont1,MyCont2,MyCont3)
writeData(workbook, "Results", MyContTime, startCol = 1, startRow = nrow(FixedOut)+nrow(RandomOut)+nrow(CorrOut)+nrow(MyCont) + 6,headerStyle = hs2)

#Testing for differences in the influence of preysize on in change in number of attackers with increasing time among the three species
MyCont1 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesdumicola:timepoint_z:preysize_z", Eff2 = "speciessarasinorum:timepoint_z:preysize_z")
MyCont2 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesmimosarum:timepoint_z:preysize_z", Eff2 = "speciessarasinorum:timepoint_z:preysize_z")
MyCont3 <- FixedEffectContrasts(model = MyModel, Eff1 = "speciesmimosarum:timepoint_z:preysize_z", Eff2 = "speciesdumicola:timepoint_z:preysize_z")
MyCont3way <- rbind(MyCont1,MyCont2,MyCont3)
writeData(workbook, "Results", MyCont3way, startCol = 1, startRow = nrow(FixedOut)+nrow(RandomOut)+nrow(CorrOut)+nrow(MyCont)+nrow(MyContTime) + 7,headerStyle = hs2)


#openXL(workbook)
#Print
saveWorkbook(workbook, file= paste("results/","Table S1 Nattackers_withnest.xlsx",sep = ""),overwrite = TRUE)
