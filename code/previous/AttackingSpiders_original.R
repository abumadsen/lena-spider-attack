#20/03/19
#Mads
#Analysis of prey size preferences - Attacking spiders

rm(list = ls())

library(ggplot2);library(tidyverse);library(lme4);library(boot);library(doBy);library(Hmisc);library(car)

#For rescaling!
#https://stats.stackexchange.com/questions/209784/rescale-predictions-of-regression-model-fitted-on-scaled-predictors
#http://stackoverflow.com/questions/10287545/backtransform-scale-for-plotting

xdata5 = read.csv("SourceData/stegopreysizeprefattackers.csv", header=T, sep=",")

#######################################
###--- PREPPING
#######################################

use <- xdata5 %>%
  select(attackers,species,id,groupsize,preysize,timepoint) %>%
  filter(!is.na(preysize)) %>%
  filter(!is.na(attackers)) %>%
  mutate_at(c("groupsize","preysize","attackers","timepoint"), funs(as.numeric(.))) %>%
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate(groupsize_z = scale(groupsize)[,1]) %>%
  mutate(timepoint_z = scale(timepoint)[,1]) %>%
  mutate(species2 = paste("S.", species)) %>%
  mutate_at(c("species","species2","id"), funs(factor(.)))

str(use)
summary(use)

#######################################
###--- INSPECTING
#######################################

hist(use$attackers)
hist(use$groupsize)
hist(use$preysize)

table(use$id,use$timepoint) #some ids have two or three trials
#I create a new id that is prey specific also
use$id_prey <- as.factor(paste(use$id, use$preysize, sep = "_"))
table(use$id_prey,use$timepoint) #some ids have two or three trials

table(use$id,use$groupsize)
table(use$species,use$groupsize)
table(use$species,use$preysize)

pdf("Final/Supplementary figure 2.pdf", family = "Times", width = 6.7, height = 5)
ggplot(use, aes(x = preysize)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(.~species2, labeller = label_bquote(col = italic(.(as.character(species2))))) +
  xlab("Preyzise (mm)") + 
  ylab("Count") +
  theme_light()
dev.off()

#Inspect change in number of attackers
ggplot(use[use$id_prey %in% levels(use$id_prey)[seq(1,20,1)],], aes(x = timepoint, y = attackers)) +
  geom_point() +
  facet_wrap(~id_prey) + 
  xlab("Preyzise (mm)") + 
  ylab("Count") +
  theme_light()

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
plot(Means$preysize_round, Means$attackers.mean, col = c(Means$species)) 	#Looks fine
#Not power to say it is not linear

#Linear effect of groupsize?
Means <- summaryBy(attackers ~ species + groupsize + timepoint, FUN = mean, data = use)
pdf("TempData/groupsize_effect.pdf")
plot(Means$groupsize, Means$attackers.mean, col = c(Means$species), cex = (Means$timepoint+1)/3) 	#Not likely to be sig
dev.off()
#There is some evidence of colinearity between groupsize and preysize here
Means <- summaryBy(preysize ~ species + groupsize, FUN = mean, data = use)
plot(Means$groupsize, Means$preysize.mean, col = c(Means$species)) 	#Not likely to be sig


# #collinearity? 
# full = glm(attackers ~ species*timepoint_log_z*preysize_z*groupsize_z, data = use,family = "poisson")
# sqrt(vif(full))
# #Also in the vif:
# #reysize_scaled:groupsize_z : 	5.349762 
# #way too high
# 
# #--- Lets try without group size
# full_nogroup = glm(attackers ~ species*timepoint_log_z*preysize_z, data = use,family = "poisson")
# sqrt(vif(full_nogroup)) #much better without group size.

#--- I think we should skip groupsize, and then make some exploration on the side of how this factor may influence the result?


#######################################
###--- ANALYSIS - without groupsize
#######################################

#looking at the model fit
full2 = glm(attackers ~ species*timepoint_z*preysize_z
            , family = "poisson", data = use)
glm.diag.plots(full2)

#Two highly influenctial
outlier_ids <- use$id_prey[use$attackers > 20]

full3 = glm(attackers ~ species*timepoint_z*preysize_z
            , family = "poisson", data = use[!use$id_prey %in% outlier_ids,])
glm.diag.plots(full3)
#much better - I should do the analysis with and without two  datapoints

#Subsetting
usePlot <- use[!use$id_prey %in% outlier_ids,]

usePlot <- usePlot %>%
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate(groupsize_z = scale(groupsize)[,1]) %>%
  mutate(timepoint_z = scale(timepoint)[,1])

#Full model
full = glmer(attackers ~ species*timepoint_z*poly(preysize_z,2,raw = TRUE)
             + (0 + timepoint_z|id_prey) + (1|id_prey)
             , family = "poisson"
             , data = use,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

# check for overdispersal
# pearson residuals
pr = residuals(full, type="pearson")
# sum of sqares
sum.dp = sum ( pr^2 )
xdf = length (residuals(full)) - length (fixef(full))
sum.dp
#[1]  199.996
xdf
#[1] 560
# test
1-pchisq(sum.dp, xdf)
#[1] 1
# null hyp is that there is no overdispersion , so this is OK
sum.dp/xdf
#[1] 0.3289408

# fine 

full_polyout = glmer(attackers ~ species*timepoint_z*preysize_z
                     + (0 + timepoint_z|id_prey) + (1|id_prey)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full, full_polyout, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# full_polyout 10 2152.5 2196.8 -1066.2   2132.5                         
# full         14 2159.0 2221.0 -1065.5   2131.0 1.4847      4     0.8293

red = glmer(attackers ~ (species+timepoint_z+preysize_z)^2
            + (0 + timepoint_z|id_prey) + (1|id_prey)
                     , family = "poisson"
                     , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red, full, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red   9 2203.0 2243.0 -1092.5   2185.0                         
# full 14 2212.1 2274.2 -1092.1   2184.1 0.9249      5     0.9684

red2 = glmer(attackers ~ species*timepoint_z + species*preysize_z
             + (0 + timepoint_z|id_prey) + (1|id_prey)
            , family = "poisson"
            , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))
anova(red2, red, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red2  8 2150.0 2185.4 -1067.0   2134.0                         
# red   9 2150.7 2190.5 -1066.3   2132.7 1.2989      1     0.2544

red3 = glmer(attackers ~ timepoint_z + species*preysize_z
             + (0 + timepoint_z|id_prey) + (1|id_prey)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))
anova(red2, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red3  7 2148.2 2179.2 -1067.1   2134.2                         
# red2  8 2150.0 2185.4 -1067.0   2134.0 0.2671      1     0.6053

red4 = glmer(attackers ~ timepoint_z + species+preysize_z
             + (0 + timepoint_z|id_prey) + (1|id_prey)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red4, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)   
# red4  6 2153.0 2179.6 -1070.5   2141.0                           
# red3  7 2148.2 2179.2 -1067.1   2134.2   6.8      1   0.009116 **
  
red5 = glmer(attackers ~ species*preysize_z
             + (0 + timepoint_z|id_prey) + (1|id_prey)
             , family = "poisson"
             , data = usePlot,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red5, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)    
# red5  6 2222.5 2249.1 -1105.2   2210.5                             
# red3  7 2148.2 2179.2 -1067.1   2134.2 76.268      1  < 2.2e-16 ***

FINAL = red3
summary(FINAL)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    0.59376    0.08423   7.049  1.8e-12 ***
#   timepoint_z                    0.29407    0.02742  10.726  < 2e-16 ***
#   speciessarasinorum             0.08020    0.11279   0.711  0.47707    
# preysize_z                     0.15795    0.09135   1.729  0.08380 .  
# speciessarasinorum:preysize_z -0.30970    0.11682  -2.651  0.00802 ** 


#####--------------------------- Looking into if there is an effect of preysize in each species:


full_dum <- glmer(attackers ~ timepoint_z + preysize_z
                  + (0 + timepoint_z|id_prey) + (1|id_prey)
      , family = "poisson"
      , data = droplevels(usePlot[usePlot$species == "dumicola",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red_dum <- glmer(attackers ~ timepoint_z
                 + (0 + timepoint_z|id_prey) + (1|id_prey)
                  , family = "poisson"
                  , data = droplevels(usePlot[usePlot$species == "dumicola",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full_dum, red_dum, test = "Chisq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# red_dum   4 970.17 984.77 -481.09   962.17                           
# full_dum  5 968.90 987.14 -479.45   958.90 3.2724      1    0.07045 .


full_sar <- glmer(attackers ~ timepoint_z + preysize_z
                  + (0 + timepoint_z|id_prey) + (1|id_prey)
                  , family = "poisson"
                  , data = droplevels(usePlot[usePlot$species == "sarasinorum",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red_sar <- glmer(attackers ~ timepoint_z
                 + (0 + timepoint_z|id_prey) + (1|id_prey)
                 , family = "poisson"
                 , data = droplevels(usePlot[usePlot$species == "sarasinorum",]),control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full_sar, red_sar, test = "Chisq")


# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)  
# red_sar   4 1186.2 1201.5 -589.10   1178.2                          
# full_sar  5 1184.3 1203.4 -587.15   1174.3 3.895      1    0.04843 *



#######################################
###--- ANALYSIS - with groupsize
#######################################

#Only use data with groupsize
use_group <- droplevels(use[!is.na(use$groupsize) & !use$id_prey %in% outlier_ids,])

use_group <- use_group %>%
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate(groupsize_z = scale(groupsize)[,1]) %>%
  mutate(timepoint_z = scale(timepoint)[,1])


full = glmer(attackers ~ species*timepoint_z*poly(preysize_z,2) + groupsize_z
             + (1 + timepoint_z|id_prey)
             , family = "poisson"
             , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

# check for overdispersal
# pearson residuals
pr = residuals(full, type="pearson")
# sum of sqares
sum.dp = sum ( pr^2 )
xdf = length (residuals(full)) - length (fixef(full))
sum.dp
#[1]  199.996
xdf
#[1] 560
# test
1-pchisq(sum.dp, xdf)
#[1] 1
# null hyp is that there is no overdispersion , so this is OK
sum.dp/xdf
#[1] 0.3289408

# fine 

full_polyout = glmer(attackers ~ species*timepoint_z*preysize_z  + groupsize_z
                     + (1 + timepoint_z|id_prey)
                     , family = "poisson"
                     , data = use,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full, full_polyout, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# full_polyout 12 2008.2 2060.6 -992.12   1984.2                         
# full         16 2014.8 2084.6 -991.41   1982.8 1.4092      4     0.8426

red = glmer(attackers ~ (species+timepoint_z+preysize_z)^2+ groupsize_z
            + (1 + timepoint_z|id_prey)
            , family = "poisson"
            , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red, full_polyout, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red          11 2006.3 2054.3 -992.18   1984.3                         
# full_polyout 12 2008.2 2060.6 -992.12   1984.2 0.1209      1     0.7281

red2 = glmer(attackers ~ species*timepoint_z + species*preysize_z + groupsize_z
             + (1 + timepoint_z|id_prey)
             , family = "poisson"
             , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))
anova(red2, red, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red2 10 2006.4 2050.0 -993.19   1986.4                         
# red  11 2006.3 2054.3 -992.18   1984.3 2.0349      1     0.1537

red3 = glmer(attackers ~ timepoint_z + species*preysize_z
             + (1 + timepoint_z|id_prey)
             , family = "poisson"
             , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))
anova(red2, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance Chisq Chi Df Pr(>Chisq)
# red3  9 2004.4 2043.7 -993.19   1986.4                        
# red2 10 2006.4 2050.0 -993.19   1986.4     0      1      0.995

red4 = glmer(attackers ~ timepoint_z + species+preysize_z + groupsize_z
             + (1 + timepoint_z|id_prey)
             , family = "poisson"
             , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red4, red3, test = "Chissq")
# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)   
# red4  8 2009.6 2044.5 -996.80   1993.6                            
# red3  9 2004.4 2043.7 -993.19   1986.4 7.2183      1   0.007216 **

red5 = glmer(attackers ~ species*preysize_z  + groupsize_z
             + (1 + timepoint_z|id_prey)
             , family = "poisson"
             , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red5, red3, test = "Chissq")
# Df    AIC    BIC   logLik deviance Chisq Chi Df Pr(>Chisq)    
# red5  8 2042.2 2077.1 -1013.11   2026.2                            
# red3  9 2004.4 2043.7  -993.19   1986.4 39.83      1  2.771e-10 ***

red6 = glmer(attackers ~ timepoint_z + species*preysize_z
             + (1 + timepoint_z|id_prey)
             , family = "poisson"
             , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))
anova(red6, red3, test = "Chissq")

# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# red6  8 2007.1 2042.0 -995.54   1991.1                           
# red3  9 2004.4 2043.7 -993.19   1986.4 4.6843      1    0.03044 *

FINAL_group = red3
summary(FINAL_group)

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)                    0.61164    0.08077   7.573 3.66e-14 ***
#   timepoint_z                    0.21827    0.03151   6.926 4.34e-12 ***
#   speciessarasinorum             0.07010    0.10744   0.652  0.51414    
# preysize_z                     0.13451    0.08816   1.526  0.12708    
# speciessarasinorum:preysize_z -0.31114    0.11330  -2.746  0.00603 ** 


#####--------------------------- Looking into if there is an effect of preysize in each species:


full_dum <- glmer(attackers ~ timepoint_z + preysize_z + groupsize_z
                  + (1 + timepoint_z|id_prey)
                  , family = "poisson"
                  , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red_dum <- glmer(attackers ~ timepoint_z  + groupsize_z
                 + (1 + timepoint_z|id_prey)
                 , family = "poisson"
                 , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full_dum, red_dum, test = "Chisq")

# Df    AIC  BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)  
# red_dum   5 2009.2 2031 -999.61   1999.2                           
# full_dum  7 2008.5 2039 -997.24   1994.5 4.7362      2    0.09366 .


full_sar <- glmer(attackers ~ timepoint_z + preysize_z  + groupsize_z
                  + (1 + timepoint_z|id_prey)
                  , family = "poisson"
                  , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red_sar <- glmer(attackers ~ timepoint_z  + groupsize_z
                 + (1 + timepoint_z|id_prey)
                 , family = "poisson"
                 , data = use_group,control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(full_sar, red_sar, test = "Chisq")


# Df    AIC    BIC  logLik deviance  Chisq Chi Df Pr(>Chisq)
# red_sar   6 2008.8 2034.9 -998.37   1996.8                         
# full_sar  7 2008.5 2039.0 -997.24   1994.5 2.2626      1     0.1325


#######################################
###--- GRAPH
#######################################

summary(FINAL)
summary(usePlot)
#Data for fitting
NEWDATdum <- expand.grid(species = "dumicola", preysize_z = seq(min(usePlot$preysize_z[usePlot$species == "dumicola"]),max(usePlot$preysize_z[usePlot$species == "dumicola"]),0.05), timepoint_z = unique(usePlot$timepoint_z))
NEWDATsar <- expand.grid(species = "sarasinorum", preysize_z = seq(min(usePlot$preysize_z[usePlot$species == "sarasinorum"]),max(usePlot$preysize_z[usePlot$species == "sarasinorum"]),0.05), timepoint_z = unique(usePlot$timepoint_z))

#All the way:
#NEWDATdum <- expand.grid(species = "dumicola", preysize_z = seq(min(use$preysize_z[use$species == "dumicola"]),max(use$preysize_z[use$species == "dumicola"]),0.05), timepoint_log = unique(use$timepoint_log))
#NEWDATsar <- expand.grid(species = "sarasinorum", preysize_z = seq(min(use$preysize_z[use$species == "sarasinorum"]),max(use$preysize_z[use$species == "sarasinorum"]),0.05), timepoint_log = unique(use$timepoint_log))


NEWDAT <- rbind(NEWDATdum, NEWDATsar)

#Predicting
NEWDAT$Fit <- predict(FINAL, NEWDAT, re.form=NA, type = "link")
head(NEWDAT)

#Rescaling preysize:
NEWDAT$preysize <- NEWDAT$preysize_z*attr(scale(use$preysize), 'scaled:scale') + attr(scale(use$preysize), 'scaled:center')

#Getting confidence interval
bb <- bootMer(FINAL, FUN = function(x)predict(x, NEWDAT, re.form=NA, type = "link"), nsim = 1000, use.u = FALSE, type = "parametric")
NEWDAT$LowC <- apply(bb$t, 2, quantile, 0.025)
NEWDAT$HighC <- apply(bb$t, 2, quantile, 0.975)

pdf("Final/Attackers_noOutliers_noGroup.pdf", family = "Times", width = 6.7, height = 6.7)

par(mfrow=c(2,2), mar = c(2.9,2.5,0.2,0.2))

mycex = 1.2


for(i in unique(NEWDAT$timepoint)){
  
  NEWDAT2 <- NEWDAT[NEWDAT$timepoint == i,]
  
  plot(NEWDAT2$preysize[NEWDAT2$species == "dumicola"], exp(NEWDAT2$Fit[NEWDAT2$species == "dumicola"]), type = "l", col = "black", ylim = c(0,10), ylab = "Number of attackers", xlab = "Preysize (mm)", cex.main = 0.7, xlim = c(0,65), las = 1, mgp = c(1.4,0.4,0), tck = -0.02,cex.lab = mycex, cex.axis = mycex, yaxt = "n")
  
  mtext(paste("Timepoint: ", round(exp(i))-1, " min", sep = ""), side = 3, line = -1.5, at = 11.5,cex = mycex-0.3)
  axis(side = 2, at = seq(0,10,2), cex.axis = mycex, las = 1, tck = -0.02, mgp = c(1.4,0.4,0))
  
  polygon(c(NEWDAT2$preysize[NEWDAT2$species == "dumicola"], rev(NEWDAT2$preysize[NEWDAT2$species == "dumicola"])), c(exp(NEWDAT2$LowC[NEWDAT2$species == "dumicola"]), rev(exp(NEWDAT2$HighC[NEWDAT2$species == "dumicola"]))), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  
  lines(NEWDAT2$preysize[NEWDAT2$species == "sarasinorum"], exp(NEWDAT2$Fit[NEWDAT2$species == "sarasinorum"]), type = "l", col = rgb(40,40,40,80, maxColorValue = 255), ylim = c(0,9))
  
  polygon(c(NEWDAT2$preysize[NEWDAT2$species == "sarasinorum"], rev(NEWDAT2$preysize[NEWDAT2$species == "sarasinorum"])), c(exp(NEWDAT2$LowC[NEWDAT2$species == "sarasinorum"]), rev(exp(NEWDAT2$HighC[NEWDAT2$species == "sarasinorum"]))), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  
  legend(-3,9.5, legend = c(expression(italic("S. dumicola")), expression(italic("S. sarasinorum"))), lty = c(1,1), col = c("black",rgb(40,40,40,80, maxColorValue = 255)), bty = "n", cex = mycex-0.1)
  legend(0,9.5, legend = c("", ""), bty = "n", cex = mycex-0.1, pch = c(16,1))
  
  #---------------		Plotting observed
  library(Hmisc)
  
  #- grouping preysize
  usePlot$preysize10 <- round(usePlot$preysize/10)
  MyMeans <- summaryBy(attackers + preysize ~ preysize10 + species, data = usePlot[usePlot$timepoint_z == i,], FUN = c(mean, sd, length))
  
  MyMeans$attackers.se <- MyMeans$attackers.sd/sqrt(MyMeans$attackers.length) 	#Laver Standard errors
  MyMeans$attackers.max <- MyMeans$attackers.mean + MyMeans$attackers.se	#Laver max for error bars
  MyMeans$attackers.min <-MyMeans$attackers.mean - MyMeans$attackers.se	#Laver min for error bars
  
  errbar(x = MyMeans$preysize.mean[MyMeans$species == "dumicola"], y = MyMeans$attackers.mean[MyMeans$species == "dumicola"], yplus = MyMeans$attackers.max[MyMeans$species == "dumicola"], yminus = MyMeans$attackers.min[MyMeans$species == "dumicola"], xlab = "", ylab = "", cap = 0.02, las = 1, mgp = c(2,0.3,0), cex.lab = mycex, cex = mycex, cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = rgb(80,80,80,80, maxColorValue = 255))
  
  errbar(x = MyMeans$preysize.mean[MyMeans$species == "sarasinorum"], y = MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], yplus = MyMeans$attackers.max[MyMeans$species == "sarasinorum"], yminus = MyMeans$attackers.min[MyMeans$species == "sarasinorum"], xlab = "", ylab = "", cap = 0.02, las = 1, mgp = c(2,0.3,0), cex.lab = mycex, cex = mycex, cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = rgb(80,80,80,80, maxColorValue = 255), pch = 1)
  
  
  #points(MyMeans$preysize.mean[MyMeans$species == "dumicola"], MyMeans$attackers.mean[MyMeans$species == "dumicola"], col = "black", cex = mycex)
  #points(MyMeans$preysize.mean[MyMeans$species == "sarasinorum"], MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], col = "red", cex = mycex)
  
  #points(use$preysize[use$species == "dumicola" & use$timepoint_log == i], use$attackers[use$species == "dumicola" & use$timepoint_log == i], pch = 19, col = rgb(50,50,50,100,maxColorValue = 255), cex = mycex-0.2)
  #points(use$preysize[use$species == "sarasinorum" & use$timepoint_log == i], use$attackers[use$species == "sarasinorum" & use$timepoint_log == i], pch = 17, col = rgb(50,50,50,100,maxColorValue = 255), cex = mycex-0.2)	
}
dev.off()




#######################################
###--- GRAPH 2
#######################################


pdf("Final/Attackers_noOutliers_noGroup_acrosstime.pdf", family = "Times", width = 6.7, height = 6.7)

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




###### WHAT ABOUT GROUP SIZE

#-- TRY TO INCLUDE OR MAYBE CREATE 4 levels AND ADD AS RANDOM



red4_gr = glmer(attackers ~ timepoint_z*groupsize_z + species*preysize_z
             + (1 + timepoint_z|id_prey)
             , family = "poisson"
             , data = use[!use$id_prey %in% outlier_ids  & !is.na(use$groupsize),],control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

red5_gr = glmer(attackers ~ timepoint_z + species*preysize_z
                + (1 + timepoint_z|id_prey)
                , family = "poisson"
                , data = use[!use$id_prey %in% outlier_ids & !is.na(use$groupsize),],control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e6)))

anova(red4_gr, red5_gr, test = "Chissq")


