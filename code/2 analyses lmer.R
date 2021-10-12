#24/06/20
#Mads F. Schou
#Analysis of prey size preferences - Attacking spiders

rm(list = ls())

pacman::p_load("ggplot2","tidyverse","doBy","lme4","boot","MCMCglmm")

use = read.csv("interm/social attack data for Mads 22.05.20_prepped.csv", header=T, sep=",")

#######################################
###--- ANALYSIS - first steps
#######################################

#Linear effect of timpoint?
Means <- summaryBy(attackers ~ species + timepoint, FUN = mean, data = use)
plot(Means$timepoint, Means$attackers.mean) 	#Maybe it needs a log transform to become linear

plot(log(Means$timepoint+1), Means$attackers.mean) 	#Yep this looks a bit better
#However the results dont differ between log and no log, so I stick to not log to ease interpretation
#use$timepoint_log <- log(use$timepoint +1)

#Linear effect of preysize?
use$preysize_round <- round(use$preysize/10)
Means <- summaryBy(attackers ~ species + preysize_round, FUN = c(mean,length), data = use)
plot(Means$preysize_round, Means$attackers.mean, col = c(Means$species)) 	#Does indicate some non-linearity - but not strong!
#Not power to say it is not linear



#looking at the model fit
full2 = glm(attackers ~ species*timepoint_z*preysize_z
            , family = "poisson", data = use)
glm.diag.plots(full2)

#Two highly influenctial
outlier_ids <- use$trial[use$attackers > 20]

full3 = glm(attackers ~ species*timepoint_z*preysize_z
            , family = "poisson", data = use[!use$trial %in% outlier_ids,])
glm.diag.plots(full3)
#much better - I should do the analysis with and without two  datapoints

#Subsetting
usePlot <- use[!use$trial %in% outlier_ids,]

usePlot <- usePlot %>%
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate(timepoint_z = scale(timepoint)[,1])



#######################################
###--- ANALYSIS 
#######################################

#Full model
full = glmer(attackers ~ species*timepoint_z*poly(preysize_z,2,raw = TRUE)
             + (0 + timepoint_z|trial) + (1|trial)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

# check for overdispersal
# pearson residuals
pr = residuals(full, type="pearson")
# sum of sqares
sum.dp = sum ( pr^2 )
xdf = length (residuals(full)) - length (fixef(full))
sum.dp
#[1]  285.0329
xdf
#[1] 748
# test
1-pchisq(sum.dp, xdf)
#[1] 1
# null hyp is that there is no overdispersion , so this is OK
sum.dp/xdf
#[1] 0.38106

# fine - maybe even overdispersion

full_polyout = glmer(attackers ~ species*timepoint_z*preysize_z
                     + (0 + timepoint_z|trial) + (1|trial)
                     , family = "poisson"
                     , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full, full_polyout, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# full_polyout 14 2572.7 2637.6 -1272.3   2544.7                         
# full         20 2579.1 2671.9 -1269.5   2539.1 5.5966      6     0.4699

red = glmer(attackers ~ (species+timepoint_z+preysize_z)^2
            + (0 + timepoint_z|trial) + (1|trial)
            , family = "poisson"
            , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red, full, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red  12 2568.9 2624.6 -1272.5   2544.9                         
# full 20 2579.1 2671.9 -1269.5   2539.1 5.8486      8     0.6642

red2 = glmer(attackers ~ species*timepoint_z + species*preysize_z
             + (0 + timepoint_z|trial) + (1|trial)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))
anova(red2, red, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red2 11 2567.8 2618.9 -1272.9   2545.8                         
# red  12 2568.9 2624.6 -1272.5   2544.9 0.9394      1     0.3324

red3 = glmer(attackers ~ timepoint_z + species*preysize_z
             + (0 + timepoint_z|trial) + (1|trial)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))
anova(red2, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red3  9 2568.3 2610.1 -1275.2   2550.3                         
# red2 11 2567.8 2618.9 -1272.9   2545.8 4.4621      2     0.1074

red4 = glmer(attackers ~ timepoint_z + species+preysize_z
             + (0 + timepoint_z|trial) + (1|trial)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red4, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# red4  7 2570.9 2603.4 -1278.5   2556.9                           
# red3  9 2568.3 2610.1 -1275.2   2550.3 6.6409      2    0.03614 *

red5 = glmer(attackers ~ species*preysize_z
             + (0 + timepoint_z|trial) + (1|trial)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red5, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# red5  8 2648.2 2685.3 -1316.1   2632.2                             
# red3  9 2568.3 2610.1 -1275.2   2550.3 81.847      1  < 2.2e-16 ***



FINAL = glmer(attackers ~ timepoint_z + species-1 + species:preysize_z
              + (0 + timepoint_z|trial) + (1|trial)
              , family = "poisson"
              , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))



summary(FINAL)

# Estimate Std. Error z value Pr(>|z|)    
# timepoint_z                    0.284601   0.025754  11.051  < 2e-16 ***
#   speciesdumicola                0.606060   0.082158   7.377 1.62e-13 ***
#   speciesmimosarum              -0.040678   0.140672  -0.289   0.7725    
# speciessarasinorum             0.732721   0.073238  10.005  < 2e-16 ***
#   speciesdumicola:preysize_z     0.148077   0.090341   1.639   0.1012    
# speciesmimosarum:preysize_z   -0.005993   0.133558  -0.045   0.9642    
# speciessarasinorum:preysize_z -0.148983   0.069556  -2.142   0.0322 *

# check for overdispersal
# pearson residuals
pr = residuals(FINAL, type="pearson")
# sum of sqares
sum.dp = sum ( pr^2 )
xdf = length (residuals(FINAL)) - length (fixef(FINAL))
sum.dp
#[1]  282.144
xdf
#[1] 748
# test
1-pchisq(sum.dp, xdf)
#[1] 1
# null hyp is that there is no overdispersion , so this is OK
sum.dp/xdf
#[1] 0.37106

#####--------------------------- Looking into if there is an effect of preysize in each species:


full_dum <- glmer(attackers ~ timepoint_z + preysize_z
                  + (0 + timepoint_z|trial) + (1|trial)
                  , family = "poisson"
                  , data = droplevels(usePlot[usePlot$species == "dumicola",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red_dum <- glmer(attackers ~ timepoint_z
                 + (0 + timepoint_z|trial) + (1|trial)
                 , family = "poisson"
                 , data = droplevels(usePlot[usePlot$species == "dumicola",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full_dum, red_dum, test = "Chisq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# red_dum   4 970.17 984.77 -481.09   962.17                           
# full_dum  5 968.90 987.14 -479.45   958.90 3.2724      1    0.07045 .


full_sar <- glmer(attackers ~ timepoint_z + preysize_z
                  + (0 + timepoint_z|trial) + (1|trial)
                  , family = "poisson"
                  , data = droplevels(usePlot[usePlot$species == "sarasinorum",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red_sar <- glmer(attackers ~ timepoint_z
                 + (0 + timepoint_z|trial) + (1|trial)
                 , family = "poisson"
                 , data = droplevels(usePlot[usePlot$species == "sarasinorum",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full_sar, red_sar, test = "Chisq")


# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)  
# red_sar   4 1186.2 1201.5 -589.10   1178.2                          
# full_sar  5 1184.3 1203.4 -587.15   1174.3 3.895      1    0.04843 *

full_mimo <- glmer(attackers ~ timepoint_z + preysize_z
                  + (0 + timepoint_z|trial) + (1|trial)
                  , family = "poisson"
                  , data = droplevels(usePlot[usePlot$species == "mimosarum",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red_mimo <- glmer(attackers ~ timepoint_z
                 + (0 + timepoint_z|trial) + (1|trial)
                 , family = "poisson"
                 , data = droplevels(usePlot[usePlot$species == "mimosarum",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full_mimo, red_mimo, test = "Chisq")

# red_mimo: attackers ~ timepoint_z + (0 + timepoint_z | trial) + (1 | trial)
# full_mimo: attackers ~ timepoint_z + preysize_z + (0 + timepoint_z | trial) + 
#   full_mimo:     (1 | trial)
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# red_mimo   4 337.14 348.55 -164.57   329.14                        
# full_mimo  5 339.13 353.39 -164.57   329.13 0.005      1     0.9436


#######################################
###--- GRAPH 2
#######################################


pdf("resukts/Attackers_noOutliers_noGroup_acrosstime.pdf", family = "Times", width = 6.7, height = 6.7)

par(mfrow=c(2,2), mar = c(2.9,2.5,0.2,0.2))

mycex = 1.2

#- grouping preysize
usePlot$preysize_cut <- as.numeric(cut_number(usePlot$preysize,4))
MyMeans <- summaryBy(preysize ~ preysize_cut + species, data = usePlot, FUN = c(mean, sd, length))
MyMeans$preysize.mean

#Rescaling timepoint:
NEWDAT$timepoint <- NEWDAT$timepoint_z*attr(scale(use$timepoint), 'scaled:scale') + attr(scale(use$timepoint), 'scaled:center')

MyPreys <- c(8,14,22,25)

for(i in c(1,2,3,4)){
  
  NEWDAT2 <- NEWDAT[floor(NEWDAT$preysize) == MyPreys[i],]
  Dum1 <- unique(NEWDAT2[NEWDAT2$species %in% "dumicola","preysize"])[1]
  Sar1 <- unique(NEWDAT2[NEWDAT2$species %in% "sarasinorum","preysize"])[1]
  NEWDAT2_sar <- NEWDAT2[NEWDAT2$species %in% "sarasinorum" & NEWDAT2$preysize == Sar1,]
  NEWDAT2_dum <- NEWDAT2[NEWDAT2$species %in% "dumicola" & NEWDAT2$preysize == Dum1,]
  NEWDAT2 <- rbind(NEWDAT2_sar,NEWDAT2_dum)
  
  plot(NEWDAT2$timepoint[NEWDAT2$species == "dumicola"], exp(NEWDAT2$Fit[NEWDAT2$species == "dumicola"]), type = "l", col = "black", ylim = c(0,5), ylab = "Number of attackers", xlab = "Timepoint", cex.main = 0.7, xlim = c(0,5), las = 1, mgp = c(1.4,0.4,0), tck = -0.02,cex.lab = mycex, cex.axis = mycex, yaxt = "n")
  
  mtext(paste("Preysize: ", MyPreys[i], " mm", sep = ""), side = 3, line = -1.5, at = 2.5,cex = mycex-0.3)
  axis(side = 2, at = seq(0,10,2), cex.axis = mycex, las = 1, tck = -0.02, mgp = c(1.4,0.4,0))
  
  polygon(c(NEWDAT2$timepoint[NEWDAT2$species == "dumicola"], rev(NEWDAT2$timepoint[NEWDAT2$species == "dumicola"])), c(exp(NEWDAT2$LowC[NEWDAT2$species == "dumicola"]), rev(exp(NEWDAT2$HighC[NEWDAT2$species == "dumicola"]))), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  
  lines(NEWDAT2$timepoint[NEWDAT2$species == "sarasinorum"], exp(NEWDAT2$Fit[NEWDAT2$species == "sarasinorum"]), type = "l", col = rgb(40,40,40,80, maxColorValue = 255), ylim = c(0,9))
  
  polygon(c(NEWDAT2$timepoint[NEWDAT2$species == "sarasinorum"], rev(NEWDAT2$timepoint[NEWDAT2$species == "sarasinorum"])), c(exp(NEWDAT2$LowC[NEWDAT2$species == "sarasinorum"]), rev(exp(NEWDAT2$HighC[NEWDAT2$species == "sarasinorum"]))), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  
  #legend(-3,1, legend = c(expression(italic("S. dumicola")), expression(italic("S. sarasinorum"))), lty = c(1,1), col = c("black",rgb(40,40,40,80, maxColorValue = 255)), bty = "n", cex = mycex-0.1)
  #legend(0,1, legend = c("", ""), bty = "n", cex = mycex-0.1, pch = c(16,1))
  
  #---------------		Plotting observed
  library(Hmisc)
  
  MyMeans <- summaryBy(attackers  ~  timepoint + species, data = usePlot[usePlot$preysize_cut == i,], FUN = c(mean, sd, length))
  
  MyMeans$attackers.se <- MyMeans$attackers.sd/sqrt(MyMeans$attackers.length) 	#Laver Standard errors
  MyMeans$attackers.max <- MyMeans$attackers.mean + MyMeans$attackers.se	#Laver max for error bars
  MyMeans$attackers.min <-MyMeans$attackers.mean - MyMeans$attackers.se	#Laver min for error bars
  
  errbar(x = MyMeans$timepoint[MyMeans$species == "dumicola"], y = MyMeans$attackers.mean[MyMeans$species == "dumicola"], yplus = MyMeans$attackers.max[MyMeans$species == "dumicola"], yminus = MyMeans$attackers.min[MyMeans$species == "dumicola"], xlab = "", ylab = "", cap = 0.02, las = 1, mgp = c(2,0.3,0), cex.lab = mycex, cex = mycex, cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = rgb(80,80,80,80, maxColorValue = 255))
  
  errbar(x = MyMeans$timepoint[MyMeans$species == "sarasinorum"], y = MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], yplus = MyMeans$attackers.max[MyMeans$species == "sarasinorum"], yminus = MyMeans$attackers.min[MyMeans$species == "sarasinorum"], xlab = "", ylab = "", cap = 0.02, las = 1, mgp = c(2,0.3,0), cex.lab = mycex, cex = mycex, cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = rgb(80,80,80,80, maxColorValue = 255), pch = 1)
  
  
  #points(MyMeans$preysize.mean[MyMeans$species == "dumicola"], MyMeans$attackers.mean[MyMeans$species == "dumicola"], col = "black", cex = mycex)
  #points(MyMeans$preysize.mean[MyMeans$species == "sarasinorum"], MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], col = "red", cex = mycex)
  
  #points(use$preysize[use$species == "dumicola" & use$timepoint_log == i], use$attackers[use$species == "dumicola" & use$timepoint_log == i], pch = 19, col = rgb(50,50,50,100,maxColorValue = 255), cex = mycex-0.2)
  #points(use$preysize[use$species == "sarasinorum" & use$timepoint_log == i], use$attackers[use$species == "sarasinorum" & use$timepoint_log == i], pch = 17, col = rgb(50,50,50,100,maxColorValue = 255), cex = mycex-0.2)	
}
dev.off()



#######################################
###--- Dealing with underdispersion
#######################################

# Underdispersion in model why? It is not overfitting but likely just the high number of 1. As there is no effect of time, we could fit max number of spiders seen within 5 min?
# Also different prey will act differently making timepoints not directly comparable?


# Note that we could leave nestID out here, as most only have one observation per nest id.

useMax <- summaryBy(attackers ~ species + nestId + preysize + trial + preysize_z + species2, data = usePlot, FUN = max)



full = glmer(attackers.max ~ species*preysize_z
                     + (1|nestId)
                     , family = "poisson"
                     , data = useMax)

red = glmer(attackers.max ~ species+preysize_z
            + (1|nestId)
            , family = "poisson"
            , data = useMax)

anova(red, full, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# red   5 880.08 896.36 -435.04   870.08                            
# full  7 872.98 895.78 -429.49   858.98 11.096      2   0.003895 **

FINAL = glmer(attackers.max ~ species-1 + species:preysize_z
      + (1|nestId)
      , family = "poisson"
      , data = useMax)


summary(FINAL)

# Estimate Std. Error z value Pr(>|z|)    
# speciesdumicola                1.09942    0.11337   9.698  < 2e-16 ***
#   speciesmimosarum               0.47996    0.17339   2.768  0.00564 ** 
#   speciessarasinorum             1.27854    0.08630  14.814  < 2e-16 ***
#   speciesdumicola:preysize_z     0.21312    0.07504   2.840  0.00451 ** 
#   speciesmimosarum:preysize_z    0.09096    0.14790   0.615  0.53852    
# speciessarasinorum:preysize_z -0.13177    0.07243  -1.819  0.06885 . 

# check for overdispersal
# pearson residuals
pr = residuals(FINAL, type="pearson")
# sum of sqares
sum.dp = sum ( pr^2 )
xdf = length (residuals(FINAL)) - length (fixef(FINAL))
sum.dp
#[1]  127.2176
xdf
#[1] 186
# test
1-pchisq(sum.dp, xdf)
#[1] 1
# null hyp is that there is no overdispersion , so this is OK
sum.dp/xdf
#[1] 0.37106
