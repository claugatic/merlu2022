
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



path.1<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/")
path.2<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/config") # config
path.3<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/input") # input
path.4<- file.path("c:/CEGM/MODEL/merluza/jjmhake/merlu2022/jjm/results") #results


# Lectura Modelos
hake1 <- jjmR::readJJM("m1", path = "config", input = "input")
