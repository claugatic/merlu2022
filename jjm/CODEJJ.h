
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


file.path()
path.1<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/")
path.2<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/config") # config
path.3<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/input") # input
path.4<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/results") #results


# Correr casos
modh1_0.00 <- jjmR::runit(mod="m1", est=T, pdf=T,
exec=file.path(path.1,"jjms"),
path=path.2,
input=path.3,
output=path.4)
