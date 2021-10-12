#24/06/20
#Mads F. Schou
#Analysis of prey size preferences - Attacking spiders

rm(list = ls())

pacman::p_load(tidyverse)

dat = read.csv("input/social attack data for Mads 22.05.20.csv", header=T, sep=",")

#######################################
###--- PREPPING
#######################################

use <- dat %>%
  mutate(trial = seq(1,nrow(.))) %>% #Several trial are made on each nest, we need to be able to differentiate bewteen these
  gather(., timepoint, attackers, Nspiders.t.0:Nspiders.t.5) %>%
  mutate(timepoint = substr(timepoint,nchar(timepoint),nchar(timepoint))) %>%
  filter(!is.na(preysize)) %>%
  filter(!is.na(attackers)) %>%
  mutate_at(c("preysize","attackers","timepoint"), funs(as.numeric(.))) %>%
  mutate(preysize = ceiling(preysize)) %>% #rounding up as many have been measured only to nearest full mm anyway
  mutate(preysize_z = scale(preysize)[,1]) %>%
  mutate(timepoint_z = scale(timepoint)[,1]) %>%
  mutate(species2 = paste("S.", species)) %>%
  mutate_at(c("species","species2","nestId","trial"), funs(factor(.))) %>%
  mutate(Comments = gsub(",",";",Comments)) #TRANSLATE "," in Comments to ";"

str(use)
summary(use)

#Missing data as NA removed, but we keep for now
use[use$nestId %in% "DumHH11b" & use$trial %in% 62,]

#Outliers?
use[use$attackers >18,] #30 and 34 -> potential outliers

write.table(use, "interm/social attack data for Mads 22.05.20_prepped.csv",sep = ",", quote = FALSE, row.names = FALSE)
