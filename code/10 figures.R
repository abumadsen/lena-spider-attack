#24/06/20
#Mads F. Schou
#Graphs of prey size preferences - Attacking spiders

pacman::p_load("dplyr","tidyr","pedantics","doBy","MCMCglmm","parallel","coda","ggplot2","Hmisc","magick")

# Fig. 1: Change in # spiders over time + preysize
# Fig. 2: Sole effect of time
# Fig. S1: Histogram of preysize distribution

#Load primary analysis
PathToModel = "code/3 N attackers MCMC/"
FileName = "N attackers no outliers.nest.RData"
PathOutput = "results/"
#Load run
load(paste(PathToModel,FileName, sep = ""))

#######################################
###--- Fig. 1: Change in # spiders over time + preysize
#######################################

#Data for fitting
NEWDATdum <- expand.grid(species = "dumicola", preysize_z = seq(min(mydat$preysize_z[mydat$species == "dumicola"]),max(mydat$preysize_z[mydat$species == "dumicola"])-0.3,0.05), timepoint_z = unique(mydat$timepoint_z))
NEWDATsar <- expand.grid(species = "sarasinorum", preysize_z = seq(min(mydat$preysize_z[mydat$species == "sarasinorum"]),max(mydat$preysize_z[mydat$species == "sarasinorum"])-0.3,0.05), timepoint_z = unique(mydat$timepoint_z))
NEWDATmim <- expand.grid(species = "mimosarum", preysize_z = seq(min(mydat$preysize_z[mydat$species == "mimosarum"]),max(mydat$preysize_z[mydat$species == "mimosarum"])-0.3,0.05), timepoint_z = unique(mydat$timepoint_z))

NEWDAT <- rbind(NEWDATdum, NEWDATsar,NEWDATmim)
NEWDAT$trial = 0
NEWDAT$nestId = 0
NEWDAT$attackers = 0

NEWDAT$preysize <- NEWDAT$preysize_z * attr(scale(mydat$preysize), 'scaled:scale') + attr(scale(mydat$preysize), 'scaled:center')
NEWDAT$timepoint <- NEWDAT$timepoint_z * attr(scale(mydat$timepoint), 'scaled:scale') + attr(scale(mydat$timepoint), 'scaled:center')


MyFit <- predict(m1.3wayPlot, newdata = NEWDAT, type = "response", interval ="confidence")
NEWDAT <- cbind(NEWDAT,MyFit)

#pdf("results/Fig1.pdf", width = 6.7, height = 6.7)
jpeg("results/Fig1.jpeg", width = 6.7, height = 6.7,quality = 500, units = "in", res = 500)

par(mfrow=c(2,2), mar = c(2.9,3,4,0.2), xpd = TRUE)

mycex = 1.2

for(i in unique(NEWDAT$timepoint)){
  
  NEWDAT2 <- NEWDAT[NEWDAT$timepoint == i,]
  
  plot(NEWDAT2$preysize[NEWDAT2$species == "dumicola"], NEWDAT2$fit[NEWDAT2$species == "dumicola"], type = "l", col = "black", ylim = c(0,10), ylab = "Number of attackers", xlab = "Preysize (mm)", cex.main = 0.7, xlim = c(0,65), las = 1, mgp = c(1.4,0.4,0), tck = -0.02,cex.lab = mycex, cex.axis = mycex, yaxt = "n")
  
  if(i == 0){
    legend(0,14, legend = c(expression(italic("S. dumicola")), expression(italic("S. sarasinorum")), expression(italic("S. mimosarum"))), lty = c(1,1,1), col = c("black",rgb(40,40,40,80, maxColorValue = 255),rgb(80,40,40,80, maxColorValue = 255)), bty = "n", cex = mycex-0.1)
    legend(0,14, legend = c("", "",""), bty = "n", cex = mycex-0.1, pch = c(16,1,2), col = c(1,1, rgb(120,80,80,180, maxColorValue = 255)))
  }
  
  mtext(paste("Timepoint: ", i, " min", sep = ""), side = 3, line = -1.5, at = 11.5,cex = mycex-0.3)
  axis(side = 2, at = seq(0,10,2), cex.axis = mycex, las = 1, tck = -0.02, mgp = c(1.4,0.4,0))
  
  polygon(c(NEWDAT2$preysize[NEWDAT2$species == "dumicola"], rev(NEWDAT2$preysize[NEWDAT2$species == "dumicola"])), c(NEWDAT2$lwr[NEWDAT2$species == "dumicola"], rev(NEWDAT2$upr[NEWDAT2$species == "dumicola"])), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  polygon(c(NEWDAT2$preysize[NEWDAT2$species == "mimosarum"], rev(NEWDAT2$preysize[NEWDAT2$species == "mimosarum"])), c(NEWDAT2$lwr[NEWDAT2$species == "mimosarum"], rev(NEWDAT2$upr[NEWDAT2$species == "mimosarum"])), border = FALSE, col = rgb(120,80,80,80, maxColorValue = 255))
  
  lines(NEWDAT2$preysize[NEWDAT2$species == "sarasinorum"], NEWDAT2$fit[NEWDAT2$species == "sarasinorum"], type = "l", col = rgb(40,40,40,80, maxColorValue = 255), ylim = c(0,9))
  lines(NEWDAT2$preysize[NEWDAT2$species == "mimosarum"], NEWDAT2$fit[NEWDAT2$species == "mimosarum"], type = "l", col = rgb(80,40,40,80, maxColorValue = 255), ylim = c(0,9))
  
  polygon(c(NEWDAT2$preysize[NEWDAT2$species == "sarasinorum"], rev(NEWDAT2$preysize[NEWDAT2$species == "sarasinorum"])), c(NEWDAT2$lwr[NEWDAT2$species == "sarasinorum"], rev(NEWDAT2$upr[NEWDAT2$species == "sarasinorum"])), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  
  
  #---------------		Plotting observed
  
  #- grouping preysize
  mydat$preysize10 <- round(mydat$preysize/10)
  MyMeans <- summaryBy(attackers + preysize ~ preysize10 + species, data = mydat[mydat$timepoint == i,], FUN = c(mean, sd, length))
  
  MyMeans$attackers.se <- MyMeans$attackers.sd/sqrt(MyMeans$attackers.length) 	#Laver Standard errors
  MyMeans$attackers.max <- MyMeans$attackers.mean + MyMeans$attackers.se	#Laver max for error bars
  MyMeans$attackers.min <-MyMeans$attackers.mean - MyMeans$attackers.se	#Laver min for error bars
  
  MyMeans$preysize.se <- MyMeans$preysize.sd/sqrt(MyMeans$preysize.length) 	#Laver Standard errors
  MyMeans$preysize.max <- MyMeans$preysize.mean + MyMeans$preysize.se	#Laver max for error bars
  MyMeans$preysize.min <-MyMeans$preysize.mean - MyMeans$preysize.se	#Laver min for error bars
  
  
  #Change number of obs with a given preysize to size factor for points
  MyMeans$SizeArea <- sqrt(sqrt(pi*MyMeans$preysize.length))
  MyMeans$SizeAreaStand <- MyMeans$SizeArea/max(MyMeans$SizeArea)
  
  #vertical bars
  errbar(x = MyMeans$preysize.mean[MyMeans$species == "dumicola"], y = MyMeans$attackers.mean[MyMeans$species == "dumicola"], yplus = MyMeans$attackers.max[MyMeans$species == "dumicola"], yminus = MyMeans$attackers.min[MyMeans$species == "dumicola"], xlab = "", ylab = "", cap = 0.0, las = 1, mgp = c(2,0.3,0), cex.lab = mycex, cex = mycex*MyMeans$SizeAreaStand[MyMeans$species == "dumicola"], cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = rgb(80,80,80,80, maxColorValue = 255))
  #horizontal bars
  arrows(x0=MyMeans$preysize.min[MyMeans$species == "dumicola"], y0=MyMeans$attackers.mean[MyMeans$species == "dumicola"], x1=MyMeans$preysize.max[MyMeans$species == "dumicola"], y1=MyMeans$attackers.mean[MyMeans$species == "dumicola"], code=3, angle=90, length=0.0, col=rgb(80,80,80,80, maxColorValue = 255), lwd=1)
  
  errbar(x = MyMeans$preysize.mean[MyMeans$species == "sarasinorum"], y = MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], yplus = MyMeans$attackers.max[MyMeans$species == "sarasinorum"], yminus = MyMeans$attackers.min[MyMeans$species == "sarasinorum"], xlab = "", ylab = "", cap = 0.0, las = 1, mgp = c(2,0.3,0), cex.lab = mycex, cex = mycex*MyMeans$SizeAreaStand[MyMeans$species == "sarasinorum"], cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = rgb(80,80,80,80, maxColorValue = 255), pch = 1)
  #horizontal bars
  arrows(x0=MyMeans$preysize.min[MyMeans$species == "sarasinorum"], y0=MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], x1=MyMeans$preysize.max[MyMeans$species == "sarasinorum"], y1=MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], code=3, angle=90, length=0.0, col=rgb(80,80,80,80, maxColorValue = 255), lwd=1)
  
  errbar(x = MyMeans$preysize.mean[MyMeans$species == "mimosarum"], y = MyMeans$attackers.mean[MyMeans$species == "mimosarum"], yplus = MyMeans$attackers.max[MyMeans$species == "mimosarum"], yminus = MyMeans$attackers.min[MyMeans$species == "mimosarum"], xlab = "", ylab = "", cap = 0.0, las = 1, mgp = c(2,0.3,0), cex.lab = mycex, cex = mycex*MyMeans$SizeAreaStand[MyMeans$species == "mimosarum"], cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, col = rgb(120,80,80,180, maxColorValue = 255), errbar.col = rgb(120,80,80,80, maxColorValue = 255), pch = 2)
  #horizontal bars
  arrows(x0=MyMeans$preysize.min[MyMeans$species == "mimosarum"], y0=MyMeans$attackers.mean[MyMeans$species == "mimosarum"], x1=MyMeans$preysize.max[MyMeans$species == "mimosarum"], y1=MyMeans$attackers.mean[MyMeans$species == "mimosarum"], code=3, angle=90, length=0.0, col=rgb(80,80,80,80, maxColorValue = 255), lwd=1)
  
  #points(MyMeans$preysize.mean[MyMeans$species == "dumicola"], MyMeans$attackers.mean[MyMeans$species == "dumicola"], col = "black", cex = mycex)
  #points(MyMeans$preysize.mean[MyMeans$species == "sarasinorum"], MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], col = "red", cex = mycex)
  
  #points(use$preysize[use$species == "dumicola" & use$timepoint_log == i], use$attackers[use$species == "dumicola" & use$timepoint_log == i], pch = 19, col = rgb(50,50,50,100,maxColorValue = 255), cex = mycex-0.2)
  #points(use$preysize[use$species == "sarasinorum" & use$timepoint_log == i], use$attackers[use$species == "sarasinorum" & use$timepoint_log == i], pch = 17, col = rgb(50,50,50,100,maxColorValue = 255), cex = mycex-0.2)	
}

dev.off()

# #change to png
# tempFig<-image_read_pdf("results/Fig1.pdf")
# image_write(tempFig, "results/Fig1.png", format = "png")
# ggsave(tempFig, file="results/Fig S1.png")


#######################################
###--- Fig. 2: Sole effect of time
#######################################

#Kunne vi evt. have en enkelt graf hvor alle data (uafhaengigt af byttestoerrelse) er poolet indenfor hver art,
#hvor #spiders er plottet over tid (0, 1, 3, 5min)? 
#Blot for tydeligt at illustrere at der ingen forskel er mellem arterne i hvordan flere angriber over tid. 
#Eller vi kunne evt. lave 3 af disse grafer, en for hver af S, M og L byttestoerrelser (e.g. S = 0-15mm; M = 15-30mm; L = 30mm+). 

#Grouping preysize by Small Medium and Large
mydat$preysize_group <- "Small (<15mm)"
mydat$preysize_group[mydat$preysize >= 15] <- "Medium (15 to 30mm)"
mydat$preysize_group[mydat$preysize > 30] <- "Large (>30mm)"
table(mydat$preysize_group)
table(mydat$preysize_group,mydat$species)



# 
# #------------ version 2
# 
# MySizes <- summaryBy(preysize_z ~ preysize_group, data = mydat, FUN = mean, keep.names = T)$preysize_z
# 
# #Data for fitting
# NEWDATdum <- expand.grid(species = "dumicola", preysize_z = MySizes, timepoint_z = unique(mydat$timepoint_z))
# NEWDATsar <- expand.grid(species = "sarasinorum", preysize_z = MySizes, timepoint_z = unique(mydat$timepoint_z))
# NEWDATmim <- expand.grid(species = "mimosarum", preysize_z = MySizes, timepoint_z = unique(mydat$timepoint_z))
# 
# NEWDAT <- rbind(NEWDATdum, NEWDATsar,NEWDATmim)
# NEWDAT$trial = 0
# NEWDAT$attackers = 0
# 
# NEWDAT$preysize <- NEWDAT$preysize_z * attr(scale(mydat$preysize), 'scaled:scale') + attr(scale(mydat$preysize), 'scaled:center')
# NEWDAT$timepoint <- NEWDAT$timepoint_z * attr(scale(mydat$timepoint), 'scaled:scale') + attr(scale(mydat$timepoint), 'scaled:center')
# 
# MyFit <- predict(m1.3wayPlot, newdata = NEWDAT, type = "response", interval ="confidence")
# NEWDAT <- cbind(NEWDAT,MyFit)
# 
# #Grouping preysize by Small Medium and Large
# NEWDAT$preysize_group <- "Small"
# NEWDAT$preysize_group[NEWDAT$preysize > 15] <- "Medium"
# NEWDAT$preysize_group[NEWDAT$preysize > 30] <- "Large"
# table(NEWDAT$preysize_group)

#----------


pdf("results/Fig2.pdf", family = "Times", width = 6.7, height = 6.7/2)

par(mfrow=c(1,3), mar = c(2.9,3,4,0.2), xpd = TRUE)

mycex = 1.2

for(i in unique(mydat$preysize_group)){
  
  MyMeans <- summaryBy(attackers ~ timepoint + species, data = mydat[mydat$preysize_group %in% i,], FUN = c(mean, sd, length))
  
  MyMeans$attackers.se <- MyMeans$attackers.sd/sqrt(MyMeans$attackers.length) 	#Laver Standard errors
  MyMeans$attackers.max <- MyMeans$attackers.mean + MyMeans$attackers.se	#Laver max for error bars
  MyMeans$attackers.min <-MyMeans$attackers.mean - MyMeans$attackers.se	#Laver min for error bars

  #Change number of obs with a given attackers to size factor for points
  MyMeans$SizeArea <- sqrt(sqrt(pi*MyMeans$attackers.length))
  MyMeans$SizeAreaStand <- MyMeans$SizeArea/max(MyMeans$SizeArea)
  
  #vertical bars
  errbar(x = MyMeans$timepoint[MyMeans$species == "dumicola"], y = MyMeans$attackers.mean[MyMeans$species == "dumicola"], yplus = MyMeans$attackers.max[MyMeans$species == "dumicola"], yminus = MyMeans$attackers.min[MyMeans$species == "dumicola"], xlab = "", ylab = "", cap = 0.0, las = 1, mgp = c(2,0.4,0), cex.lab = mycex, cex.axis = mycex, xaxt = "n", tck = -0.02, errbar.col = "black", ylim = c(0,6), xlim = c(-0.3,5.3), cex =0.01, col = "white")
  points(x = MyMeans$timepoint[MyMeans$species == "dumicola"], y = MyMeans$attackers.mean[MyMeans$species == "dumicola"],bg = rgb(150,150,190,255, maxColorValue = 255), pch = 22, col = "black", cex = mycex+0.5)
  
  errbar(x = MyMeans$timepoint[MyMeans$species == "sarasinorum"]+0.1, y = MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], yplus = MyMeans$attackers.max[MyMeans$species == "sarasinorum"], yminus = MyMeans$attackers.min[MyMeans$species == "sarasinorum"], xlab = "", ylab = "", cap = 0.0, las = 1, mgp = c(2,0.4,0), cex.lab = mycex, cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = "black", cex = 0.01,col = "white")
  points(x = MyMeans$timepoint[MyMeans$species == "sarasinorum"]+0.1, y = MyMeans$attackers.mean[MyMeans$species == "sarasinorum"], pch = 23,bg = rgb(150,190,150,255, maxColorValue = 255), col = "black", cex = mycex+0.5)
  
  errbar(x = MyMeans$timepoint[MyMeans$species == "mimosarum"]-0.1, y = MyMeans$attackers.mean[MyMeans$species == "mimosarum"], yplus = MyMeans$attackers.max[MyMeans$species == "mimosarum"], yminus = MyMeans$attackers.min[MyMeans$species == "mimosarum"], xlab = "", ylab = "", cap = 0.0, las = 1, mgp = c(2,0.4,0), cex.lab = mycex, cex.axis = mycex, xaxt = "n", tck = -0.02, add = TRUE, errbar.col = "black", cex = 0.01, col = "white")
  points(x = MyMeans$timepoint[MyMeans$species == "mimosarum"]-0.1, y = MyMeans$attackers.mean[MyMeans$species == "mimosarum"], pch = 21, col = "black",bg = rgb(190,150,150,255, maxColorValue = 255), cex = mycex+0.5)
  
  if(i == "Small (<15mm)"){
    legend(0,7.67, legend = c(expression(italic("S. dumicola")), expression(italic("S. sarasinorum")), expression(italic("S. mimosarum"))), lty = c(1,1,1), col = c("black","black","black"), pt.bg = c(rgb(150,150,190,255, maxColorValue = 255),rgb(150,190,150,255, maxColorValue = 255),rgb(190,150,150,255, maxColorValue = 255)), bty = "n", cex = mycex-0.1, pch = c(22,23,21))
    #legend(0,7.6, legend = c("", "",""), bty = "n", cex = mycex-0.1, pch = c(16,1,2), col = c(1,1, rgb(120,80,80,180, maxColorValue = 255)))
  }
  
  #Format
  mtext(paste("Preysize: ", i, sep = ""), side = 3, line = -1.5, at = 2.5,cex = mycex-0.5)
  axis(side = 1, at = seq(0,5,1), cex.axis = mycex, las = 1, tck = -0.02, mgp = c(1.4,0.4,0))
  mtext("Time (min)",side = 1, cex = mycex-0.3, line = 1.6)
  mtext("Number of attackers",side = 2, cex = mycex-0.3, line = 1.5)
  
  
  
  # NEWDAT2 <- NEWDAT[NEWDAT$preysize_group == i,]
  # 
  # lines(NEWDAT2$timepoint[NEWDAT2$species == "dumicola"], NEWDAT2$fit[NEWDAT2$species == "dumicola"], type = "l", col = rgb(40,40,40,80, maxColorValue = 255), ylim = c(0,9))
  # polygon(c(NEWDAT2$timepoint[NEWDAT2$species == "dumicola"], rev(NEWDAT2$timepoint[NEWDAT2$species == "dumicola"])), c(NEWDAT2$lwr[NEWDAT2$species == "dumicola"], rev(NEWDAT2$upr[NEWDAT2$species == "dumicola"])), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  # 
  # lines(NEWDAT2$timepoint[NEWDAT2$species == "sarasinorum"], NEWDAT2$fit[NEWDAT2$species == "sarasinorum"], type = "l", col = rgb(40,40,40,80, maxColorValue = 255), ylim = c(0,9))
  # polygon(c(NEWDAT2$timepoint[NEWDAT2$species == "sarasinorum"], rev(NEWDAT2$timepoint[NEWDAT2$species == "sarasinorum"])), c(NEWDAT2$lwr[NEWDAT2$species == "sarasinorum"], rev(NEWDAT2$upr[NEWDAT2$species == "sarasinorum"])), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  # 
  # lines(NEWDAT2$timepoint[NEWDAT2$species == "mimosarum"], NEWDAT2$fit[NEWDAT2$species == "mimosarum"], type = "l", col = rgb(40,40,40,80, maxColorValue = 255), ylim = c(0,9))
  # polygon(c(NEWDAT2$timepoint[NEWDAT2$species == "mimosarum"], rev(NEWDAT2$timepoint[NEWDAT2$species == "mimosarum"])), c(NEWDAT2$lwr[NEWDAT2$species == "mimosarum"], rev(NEWDAT2$upr[NEWDAT2$species == "mimosarum"])), border = FALSE, col = rgb(80,80,80,80, maxColorValue = 255))
  # 

}

dev.off()



#######################################
###--- Fig. S1: Histogram of preysize distribution
#######################################

#png("results/Fig S1.png", family = "Times", width = 6.7, height = 5)

mydat2 <- unique(mydat[,c("trial","species2","preysize")])
mydat2$species3 <- as.character(mydat2$species2)

FigS1 <- ggplot(mydat2, aes(x = preysize)) +
  geom_histogram(binwidth = 1) +
  facet_wrap(.~species3, labeller = label_bquote(col = italic(.(species3)))) +
  xlab("Prey size (mm)") + 
  ylab("Number of trials") +
  theme_light()
#dev.off()

ggsave(FigS1, file="results/Fig S1.png" ,width = 6.7, height = 5)




#######################################
###--- Fig. 3: Time to attack
#######################################

#Load primary analysis
PathToModel = "code/5 Time to attack/"
FileName = "Time to attack.RData"
PathOutput = "results/"
#Load run
load(paste(PathToModel,FileName, sep = ""))



#Data for fitting
NEWDATdum <- expand.grid(species = "dumicola", preysize_z = seq(min(mydat$preysize_z[mydat$species == "dumicola"]),max(mydat$preysize_z[mydat$species == "dumicola"])-0.3,0.05))
NEWDATsar <- expand.grid(species = "sarasinorum", preysize_z = seq(min(mydat$preysize_z[mydat$species == "sarasinorum"]),max(mydat$preysize_z[mydat$species == "sarasinorum"])-0.3,0.05))
NEWDATmim <- expand.grid(species = "mimosarum", preysize_z = seq(min(mydat$preysize_z[mydat$species == "mimosarum"]),max(mydat$preysize_z[mydat$species == "mimosarum"])-0.3,0.05))

NEWDAT <- rbind(NEWDATdum, NEWDATsar,NEWDATmim)
NEWDAT$nestId = 0
NEWDAT$LatencyToAttackSec = 0

NEWDAT$preysize <- NEWDAT$preysize_z * attr(scale(mydat$preysize), 'scaled:scale') + attr(scale(mydat$preysize), 'scaled:center')

MyFit <- predict(m1_3.1way.nestid[[1]], newdata = NEWDAT, type = "response", interval ="confidence")
NEWDAT <- cbind(NEWDAT,MyFit)

#pdf("results/Fig1.pdf", width = 6.7, height = 6.7)
jpeg("results/Fig3.jpeg", width = 6.7, height = 5,quality = 500, units = "in", res = 500)

par(mfrow=c(1,3), mar = c(2.9,3,0.2,0.0), xpd = TRUE)

mycex = 1.2

dumiCol <- rgb(150,150,190,125, maxColorValue = 255)
saraCol <- rgb(150,190,150,125, maxColorValue = 255)
mimoCol <- rgb(190,150,150,125, maxColorValue = 255)

plot(NEWDAT$preysize[NEWDAT$species == "dumicola"], NEWDAT$fit[NEWDAT$species == "dumicola"], type = "l", lwd = 1.5, col = "white", ylim = c(0,1200), ylab = "Time to attack (min)", xlab = "Preysize (mm)", cex.main = 0.7, xlim = c(0,65), las = 1, mgp = c(1.6,0.5,0), tck = -0.02,cex.lab = mycex, cex.axis = mycex, yaxt = "n")
polygon(c(NEWDAT$preysize[NEWDAT$species == "dumicola"], rev(NEWDAT$preysize[NEWDAT$species == "dumicola"])), c(NEWDAT$lwr[NEWDAT$species == "dumicola"], rev(NEWDAT$upr[NEWDAT$species == "dumicola"])), border = FALSE, col = dumiCol)
lines(NEWDAT$preysize[NEWDAT$species == "dumicola"], NEWDAT$fit[NEWDAT$species == "dumicola"], type = "l", col = "black")
points(x = jitter(mydat$preysize[mydat$species == "dumicola"]), y = mydat$LatencyToAttackSec[mydat$species == "dumicola"], pch = 22, col = "black",bg = dumiCol, cex = mycex-0.2)
axis(side = 2, at = seq(0,1200,120), labels = seq(0,1200,120)/60, cex.axis = mycex, las = 1, tck = -0.02, mgp = c(1.4,0.4,0))
mtext(expression(italic("S. dumicola")), side = 3, line = -2)


plot(NEWDAT$preysize[NEWDAT$species == "mimosarum"], NEWDAT$fit[NEWDAT$species == "mimosarum"], type = "l", lwd = 1.5, col = "white", ylim = c(0,1200), ylab = "Time to attack (min)", xlab = "Preysize (mm)", cex.main = 0.7, xlim = c(0,65), las = 1, mgp = c(1.6,0.5,0), tck = -0.02,cex.lab = mycex, cex.axis = mycex, yaxt = "n")
polygon(c(NEWDAT$preysize[NEWDAT$species == "mimosarum"], rev(NEWDAT$preysize[NEWDAT$species == "mimosarum"])), c(NEWDAT$lwr[NEWDAT$species == "mimosarum"], rev(NEWDAT$upr[NEWDAT$species == "mimosarum"])), border = FALSE, col =mimoCol)
lines(NEWDAT$preysize[NEWDAT$species == "mimosarum"], NEWDAT$fit[NEWDAT$species == "mimosarum"], type = "l", col = "black")
points(x = jitter(mydat$preysize[mydat$species == "mimosarum"]), y = mydat$LatencyToAttackSec[mydat$species == "mimosarum"], pch = 21, col = "black",bg = mimoCol, cex = mycex-0.2)
axis(side = 2, at = seq(0,1200,120), labels = seq(0,1200,120)/60, cex.axis = mycex, las = 1, tck = -0.02, mgp = c(1.4,0.4,0))
mtext(expression(italic("S. mimosarum")), side = 3, line = -2)

plot(NEWDAT$preysize[NEWDAT$species == "sarasinorum"], NEWDAT$fit[NEWDAT$species == "sarasinorum"], type = "l", lwd = 1.5, col = "white", ylim = c(0,1200), ylab = "Time to attack (min)", xlab = "Preysize (mm)", cex.main = 0.7, xlim = c(0,65), las = 1, mgp = c(1.6,0.5,0), tck = -0.02,cex.lab = mycex, cex.axis = mycex, yaxt = "n")
polygon(c(NEWDAT$preysize[NEWDAT$species == "sarasinorum"], rev(NEWDAT$preysize[NEWDAT$species == "sarasinorum"])), c(NEWDAT$lwr[NEWDAT$species == "sarasinorum"], rev(NEWDAT$upr[NEWDAT$species == "sarasinorum"])), border = FALSE, col = saraCol)
lines(NEWDAT$preysize[NEWDAT$species == "sarasinorum"], NEWDAT$fit[NEWDAT$species == "sarasinorum"], type = "l", col = "black")
points(x = jitter(mydat$preysize[mydat$species == "sarasinorum"]), y = mydat$LatencyToAttackSec[mydat$species == "sarasinorum"], pch = 23, col = "black",bg = saraCol, cex = mycex-0.2)
axis(side = 2, at = seq(0,1200,120), labels = seq(0,1200,120)/60, cex.axis = mycex, las = 1, tck = -0.02, mgp = c(1.4,0.4,0))
mtext(expression(italic("S. sarasinorum")), side = 3, line = -2)

dev.off()


# ggplot(data = mydat, aes(x = preysize, y = LatencyToAttackSec/60), group = species) +
#   stat_smooth(aes(colour = species), method = "lm") +
#   geom_point(aes(colour = species))
