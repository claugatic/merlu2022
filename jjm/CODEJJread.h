
rm(list=ls()) # LIMPIA DATOS
WD<- getwd()  # UBICA DIRECTOR
setwd(WD)     # RECONOCE DIRECTORIO


#install.packages("devtools")
#devtools::install_github("SPRFMO/jjmR")
library(jjmR)
library(here)
library(dplyr)
library(foreach)
library(doParallel)
registerDoParallel(6)


ruta= "c:/CEGM/MODEL/J22"
file.path()
path.1<- file.path("c:/CEGM/MODEL/DJJM/jjm/src")
path.2<- file.path("c:/CEGM/MODEL/DJJM/jjm/assessment/config") # config
path.3<- file.path("c:/CEGM/MODEL/DJJM/jjm/assessment/input") # input
path.4<- file.path("c:/CEGM/MODEL/DJJM/jjm/assessment/results") #results

# Primera parte codigo
casos <- list.files(path=path.2,"*.ctl")
casos
casos <- as.character(strsplit(casos, split=".ctl"))
casosr <- casos[c(1,2,3,4,5,6,7,8,9,10,11,12,13)];
casosr



# Lectura

# compare models
# h1_1.02 h1_1.07
mods2compare <- casos[c(1,2,3)]
mods2compare
comp.plots <- compareModels(mods2compare)
plot(comp.plots,combine=T,what="biomass",stack=F,main="Biomass")
plot(comp.plots,combine=T,what="recruitment",stack=F,main="Recruitment")
plot(comp.plots,combine=T,what="ftot",stack=F,main="Fishing mortality")

# Modelos
modh1_0.00 <- jjmR::readJJM("h1_0.00", path = "config", input = "input")
modh1_0.03 <- jjmR::readJJM("h1_0.03", path = "config", input = "input")
modh1_0.04 <- jjmR::readJJM("h1_0.04", path = "config", input = "input")

# like
LL <- cbind(summary(modh1_0.00)$like,
summary(modh1_0.03)$like,
summary(modh1_0.04)$like)
LL

# diagnostic
h2.diag <- diagnostics(modh1_0.00, plot=F)
b <-plot(h2.diag, var = "fishedUnfishedBiomass")
plot(h2.diag, var = "summarySheet")
jjmR::kobe(modh1_0.00)

