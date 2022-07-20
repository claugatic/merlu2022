library(jjmR)
library(ggplot2)
library(tidyverse)
# mod0.00 <- runit(geth("0.00"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
getwd()
geth <- function(mod,h=hyp) paste0(h,"_", mod)
hyp="h1"
setwd("..")
mod0.00 <- readJJM("h1_0.00", path = "config", input = "input")
mod1.00 <- readJJM("h1_0.04", path = "config", input = "input")
mod1.03 <- readJJM("h1_1.03", path = "config", input = "input")
mod1.04 <- readJJM("h1_1.04", path = "config", input = "input")
mod1.07 <- readJJM("h1_1.07", path = "config", input = "input")
mod1.08 <- readJJM("h1_1.08", path = "config", input = "input")
mod1.09 <- readJJM("h1_1.09", path = "config", input = "input")
mod1.10 <- readJJM("h1_1.10", path = "config", input = "input")
mod1.13ls <- readJJM("h1_1.13ls", path = "config", input = "input")
mod1.13hs <- readJJM("h1_1.13hs", path = "config", input = "input")
mod1.13ll <- readJJM("h1_1.13ll", path = "config", input = "input")

#mod1.04 <- runit(geth("1.04"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
#mod1.00 <- runit(geth("0.00"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
#mod1.07 <- runit(geth("1.07"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
#mod1.08 <- runit(geth("1.08"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
#mod1.09 <- runit(geth("1.09"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
#mod1.10 <- runit(geth("1.10"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
hyp="h1"
mod0.00 <- runit(geth("0.00"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod1.11 <- runit(geth("1.11"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod1.12 <- runit(geth("1.12"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod1.13ls <- mod1.13
mod1.13ls <- runit(geth("1.13ls"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod1.13hs <- runit(geth("1.13hs"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod1.13hl<- runit(geth("1.13hl"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod1.13ll<- runit(geth("1.13ll"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")

hyp="h2"
mod2.14 <- runit(geth("1.14"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")


# mod_prev <- readJJM(geth("1.00"), path = "config", input = "input")
# save(mod_prev, file="results/mod_prev_h1.Rdat")

#load("results/mod_prev_h1.Rdat")


Mods <- combineModels(mod1.07,mod1.08)
Mods <- combineModels(mod1.07,mod1.08)
Mods <- combineModels(mod1.00,mod1.07)
Mods <- combineModels(mod1.00,mod1.04,mod1.07,mod1.08)
Mods[[1]][[5]]$Stock_1$FW_N_Chile
Mods[[2]][[5]]$Stock_1$FW_N_Chile
Mods[[1]][[5]]$Stock_1$FW_SC_Chile_PS
Mods[[2]][[5]]$Stock_1$FW_SC_Chile_PS
Mods <- combineModels(mod1.00,mod1.10,mod1.11,mod1.12,mod1.13)
summary(Mods)
plot(Mods,combine=TRUE,what='ssb',stack=F)
plot(Mods,combine=TRUE,what='recruitment',stack=F)
library(tidyverse)
df1<-as_tibble(mod1.13ls[[1]][[5]][[1]]$df1)
df2<-as_tibble(mod1.13hs[[1]][[5]][[1]]$df1)
df3<-as_tibble(mod1.13ll[[1]][[5]][[1]]$df1)
df4<-as_tibble(mod1.13hl[[1]][[5]][[1]]$df1)

names(df1)<-names(df2)<-c("type","year","est","std","lb","ub")
df<-rbind(df1 |> mutate(model="Low"),
          df2 |> mutate(model="High"))
names(df4)<-names(df3)<-names(df1)<-names(df2)<-c("type","year","est","std","lb","ub")
df<-rbind(df1 |> mutate(model="Low, short"),
          df2 |> mutate(model="High, short"),
          df3 |> mutate(model="Low, long"),
          df4 |> mutate(model="High, long"))

df <- df %>% mutate(
  year=as.numeric(year),
  est=as.numeric(est),
  std=as.numeric(std),
  lb =as.numeric(lb),
  ub =as.numeric(ub)
)
unique(df$type)
blim <- df |> filter(type=="SSB0_Dynamic") |>  group_by(model) |> summarise(min(est))
blim
df %>% 
  ggplot(aes(x=year,y=est,ymin=lb,ymax=ub,color=model,fill=model)) + 
  geom_ribbon(alpha=.4)+ geom_line(size=.5,alpha=.5) + expand_limits(y=0) + ggthemes::theme_few() + ylab("Estimates") +
  facet_wrap(type~.,scales="free")

#df %>% filter(type=="Total_Biom") %>%

df %>%
  ggplot(aes(x=year,y=est,ymin=lb,ymax=ub)) + geom_ribbon(fill="salmon") + expand_limits(y=0) + theme_bw() + ylab("Estimates") +
  facet_wrap(type~.,scales="free")

#df |> filter(type=="SSB0_Dynamic") |> 
df |> filter(type=="SSB_Nofishing") #|> 
  ggplot(aes(x=year,y=est,ymax=ub,ymin=lb)) + geom_ribbon()

selectivities <- get_selectivities(Mods)
compareModels(Mods)

#Renzo's
mod=mod1.10$h1_1.10

rm(weights.table, Ws)
W_labels <- names(mod[[1]][c(3:8)])
##
for(i in 1:length(W_labels)){
  ##
  if(i==1)
    { Ws <- as.numeric(mod[[1]][c(i+2)])}
  else 
    {Ws <- rbind(Ws, as.numeric(mod[[1]][c(i+2)]))}
  ##
}
##
weights.table <- Ws
rownames(weights.table)  <- W_labels; colnames(weights.table) <- 'FWs'
weights.table
rownames(weights.table)  <- W_labels; colnames(weights.table) <- 'FWs'
weights.table




plot_selectivities(selectivities)

plot(Mods,combine=TRUE,what="selectivity",stack=F,fleet="fsh")

plot(Mods,what="Recruitment")

mod_diag <- diagnostics(mod1.03,plot=F)

mod <- mod1.08


rm(weights.table, Ws)
fw_table <-function(mod=mod1.08) {
  W_labels <- names(mod[[1]][c(3:8)])
  for(i in 1:length(W_labels)){
    if(i==1){ Ws <- as.numeric(mod[[1]][c(i+2)])}
    else {Ws <- rbind(Ws, as.numeric(mod[[1]][c(i+2)]))}
  }
rownames(weights.table)  <- W_labels; colnames(weights.table) <- 'FWs'
return(weights.table)
}
fw_table()
plot(mod_diag,var="catchResidualsByFleet")
plot(mod_diag,var="absoluteResidualCatchByFleet")
plot(diagnostics(mod1.03,plot=F),var="residualsCatchAtAgeByFleet")
plot(diagnostics(mod1.08,plot=F),var="residualsCatchAtAgeByFleet")
plot(mod_diag,var="residualsCatchAtLengthByFleet")
plot(diagnostics(mod1.07,plot=F),var="ageFitsCatch")
plot(diagnostics(mod1.07,plot=F),var="fisheryMeanAge")
plot(diagnostics(mod1.08,plot=F),var="fisheryMeanAge")
plot(diagnostics(mod1.08,plot=F),var="ageFitsCatch")
plot(diagnostics(mod1.08,plot=F),var="ageFitsCatch$SC_Chile_PS")

plot(diagnostics(mod1.07,plot=F),var="predictedObservedIndices")
plot(diagnostics(mod1.08,plot=F),var="predictedObservedIndices")




  library(jjmR)
  library(ggthemes)
  library(ggridges)
  library(tidyverse)
  library(PBSadmb)
  source("../R/read.admb.R")
  source("../R/read-admb.R")
p2<-readList("newAge/For_R_1.rep")

library(jjmR)

          mods2compare <- c("0.03","0.04")
          obj <- compareModels(geth(mods2compare,"h2"))
          plot(obj,combine=T,what="biomass",stack=F,main="Biomass")
          
hyp <- "h2"
mod0.03 <- runit(geth("0.03"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod0.04 <- runit(geth("0.04"),pdf=TRUE,portrait=F,est=TRUE,exec="../src/jjms")
mod0.03 <- readJJM(geth("0.03"),,exec="../src/jjms")
mod0.04 <- readJJM(geth("0.04"),pdf=TRUE,portrait=F,exec="../src/jjms")
h4.mod <- readJJM(h1nm, path = "config", input = "input")
h3.diag <- diagnostics(mod0.03,plot=F)
h4.diag <- diagnostics(mod0.04,plot=F)
h4.diag


p4<-readList("sel4/For_R_1.rep")
p2<-readList("For_R_2.rep")

p1<-readList("For_R_1.rep")
sf1<-as_tibble(rbind(p1$sel_fsh_1,
p1$sel_fsh_2,
p1$sel_fsh_3,
p1$sel_fsh_4))

sf2<-as_tibble(rbind(p2$sel_fsh_1,
p2$sel_fsh_2,
p2$sel_fsh_3,
p2$sel_fsh_4))

sf2<-as_tibble(rbind(p4$sel_fsh_1,
sf2<-as_tibble(rbind(p3$sel_fsh_1,
p3$sel_fsh_2,
p3$sel_fsh_3,
p3$sel_fsh_4))
sf2<-as_tibble(rbind(p4$sel_fsh_1,
p4$sel_fsh_2,
p4$sel_fsh_3,
p4$sel_fsh_4))
names(sf1) <- c("idx","year",1:12)
names(sf2) <- c("idx","year",0:12)
bind_rows(sf1,sf2)%>% arrange_col()
sf1 %>% pivot_longer() 
selfsh <- rbind(sf1%>%mutate(model="base"),sf2%>%mutate(model="newAge")) %>% mutate(idx=ifelse(idx==4,3,idx))

selfsh <- rbind(sf1%>%mutate(model="base"),sf2%>%mutate(model="3param_sel")) %>% mutate(idx=ifelse(idx==4,3,idx))

selfsh <- rbind(sf1%>%mutate(model="base"),sf2%>%mutate(model="3param_sel")) %>% mutate(idx=ifelse(idx==4,3,idx))
selfsh$idx <- p1$Fshry_names[selfsh$idx]

selfsh %>%
  gather(age,selectivity,3:14)  %>% filter(year>=1980) %>% mutate(age=as.numeric(age),idx=as.factor(idx))  %>%
  ggplot(aes(x=age,y=fct_rev(as.factor(year)),height = selectivity,fill=model,color=model,alpha=.3)) +
     geom_density_ridges(stat = "identity",scale=3,alpha = .3) + ylab("Year")+ theme_few() + 
     theme(legend.position="none")+facet_grid(model~idx) 

######################################
## ADNUTS try
library(adnuts)
# Model name
m <- './jjms'
# Directory  
d <- 'test'
# Assumes a converged MLE model has already been run....
setwd(d)
system(paste(m, '-nox -binp jjms.bar -phase 22 -mcmc 10 -hbf 1'))
setwd('..')

## Two different ways to gets NUTS working. First is to use the
## Hessian (metric) just like with the RMW. Note the control argument.
## Two different ways to gets NUTS working. First is to use the
## Hessian (metric) just like with the RMW. Note the control argument.
iter <- 2000 # maybe too many...depends are number cores...I used 8...
chains=8
fit.mle <- sample_nuts(model=m, path=d, iter=iter, warmup=iter/4, 
                   chains=chains, cores=chains, control=list(metric='mle'))
summary(fit.mle)
#saveRDS(fit.mle,"fit.mle.RDS")
# now test for 3-parameter logistic
d <- 'test/sel3'
iter <- 1000 # maybe too many...depends are number cores...I used 8...
chains<-8

fit.3par <- sample_nuts(model=m, path=d, iter=iter, warmup=iter/4,
                   chains=chains, cores=chains, control=list(metric='mle'))
saveRDS(fit.3par,"fit.3par.RDS")
d <- 'test/sel4'
fit.3par.sel4 <- sample_nuts(model=m, path=d, iter=iter, warmup=iter/4,
                   chains=chains, cores=chains, control=list(metric='mle'))
######################################
mcldf <- read_csv("mclike.csv")
mcldf
mnssb <- mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% group_by(type) %>% summarise(mean_SSB=mean(SSB)) 
mnssb
bin <- seq(5000,12000,by=500); bin
mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% mutate(bin=cut(SSB, breaks=bin, labels=FALSE)) %>% group_by(type) %>% mutate(NLL=value-mean(value)) %>%
ungroup() %>% group_by(type,bin) %>% summarise(NLL=mean(NLL))%>%
ggplot(aes(x=bin,y=NLL,color=type)) + theme_few()+ geom_line() #+ facet_wrap(type~.) 

mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% group_by(type) %>% mutate(NLL=value-mean(value)) %>%
ggplot(aes(x=SSB,y=NLL,color=type)) + theme_few()+ geom_point(alpha=.2) + facet_wrap(type~.,scales="free") + stat_smooth()

#mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% mutate(bin=cut(SSB, breaks=bin, labels=FALSE)) %>% group_by(type) %>% mutate(NLL=value-min(value)) %>%
tmp <- mcldf %>% filter(stock==1,type!="ind_len",type!="priors") #%>% mutate(SSB=arules::discretize(SSB, breaks = 3, labels = c("Low","Medium","High"))
tmp$SSB <- arules::discretize(tmp$SSB, breaks = 5, labels = c("Low","Med-low","Medium","Med-high","High"))
tmp %>% group_by(type) %>% mutate(NLL=value-min(value)) %>% group_by(type,SSB) %>% summarise(NLL=mean(NLL)) %>%
ggplot(aes(x=SSB,y=NLL,color=type)) + theme_few()+ geom_point(size=2) + facet_wrap(type~.) 
# mutate(bin=cut(SSB, breaks=bin, labels=FALSE)) 

pairs_admb(fit.mle, pars=1:6, order='slow')
pairs_admb(fit.3par, pars=1:6, order='slow')
pairs_admb(fit.3par.sel4, pars=1:6, order='slow')
pairs_admb(fit.mle, pars=1:6, order='fast')
print(fit.mle)
plot_sampler_params(fit.mle)
launch_shinyadmb(fit.mle)
launch_shinyadmb(fit.3par)
library(tidyverse)
mcdf <- read.table("test/mceval.rep",header=TRUE)
mcdf <- as.tibble(mcdf)
glimpse(mcdf)
unique(mcdf$type)
unique(mcdf$Age)
mcdf %>% filter(Year>1990,type=="SSB", Age=="all_stock_1") %>% mutate(Year=as.factor(Year) ) %>%
  ggplot(aes(x=Year,y=value)) + geom_violin(color="salmon",fill="salmon") + ggthemes::theme_few() 

mcdf %>% filter(Year>1990,type=="Recruits", Age=="Age_1_stock_1") %>% mutate(Year=as.factor(Year) ) %>%
  ggplot(aes(x=Year,y=value)) + geom_violin(color="salmon",fill="salmon") + ggthemes::theme_few() 

mcdf %>% filter(Year>1990,type=="Depletion", Age=="all_stock_1") %>% mutate(Year=as.factor(Year) ) %>%
  ggplot(aes(x=Year,y=value)) + geom_violin(color="salmon",fill="salmon") + ggthemes::theme_few() 
  
 # + xlim(c(0,8000))

mon <- monitor(fit.mle$samples, warmup=fit.mle$warmup, print=FALSE)


r

df <-rbind( data.frame(
    Year=p1$SSB[,1],
    SSB=p1$SSB[,2],
    lb=p1$SSB[,4],
    ub=p1$SSB[,5],
    model="Peru (Dec)"),
  data.frame(
    Year=p2$SSB[,1],
    SSB=p2$SSB[,2],
    lb=p2$SSB[,4],
    ub=p2$SSB[,5],
    model="SPRFMO SC08")
  )
df <-rbind( data.frame(
    Year=p1$R[,1],
    Recruits=p1$R[,2],
    lb=p1$R[,4],
    ub=p1$R[,5],
    model="Peru (Dec)"),
  data.frame(
    Year=p2$R[,1],
    Recruits=p2$R[,2],
    lb=p2$R[,4],
    ub=p2$R[,5],
    model="SPRFMO SC08")
  )
glimpse(df)
df %>% ggplot(aes(x=Year,y=SSB,ymin=lb,ymax=ub,color=model,fill=model)) + geom_ribbon(alpha=.3) +
geom_line(size=2)+ theme_few()
df %>% ggplot(aes(x=Year,y=Recruits,ymin=lb,ymax=ub,color=model,fill=model)) + geom_ribbon(alpha=.3) +
geom_line(size=2)+ theme_few()
read.admb <- function(ifile) { 
  ret=read.fit(ifile)
  
  fn=paste(ifile,'.rep', sep='')
  A=read.rep(fn)
  A$fit=ret
  
  pfn=paste(ifile,'.psv',sep='')
  if(file.exists(pfn))
    A$post.samp=read.psv(pfn)
  
  return(A)
}

read.fit <- function(ifile) {
  # __Example:             
  # file <-("~/admb/simple")
  # A <- reptoRlist(file)
  # Note there is no extension on the file name.
  
  ## The following is a contribution from:
  ## Anders Nielsen that reads the par & cor files.
  ret<-list() 
  parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
   what='', n=16, quiet=TRUE)[c(6,11,16)]) 
  ret$nopar<-as.integer(parfile[1]) 
  ret$nlogl<-parfile[2] 
  ret$maxgrad<-parfile[3] 
  file<-paste(ifile,'.cor', sep='') 
  lin<-readLines(file) 
  ret$npar<-length(lin)-2 
  ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
  sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
  ret$names<-unlist(lapply(sublin,function(x)x[2])) 
  ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
  ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
  ret$cor<-matrix(NA, ret$npar, ret$npar) 
  corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
  ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
  ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
  ret$cov<-ret$cor*(ret$std%o%ret$std)
  return(ret)
}

read.rep <- function(fn) {
  # The following reads a report file
  # Then the 'A' object contains a list structure
  # with all the elemements in the report file.
  # In the REPORT_SECTION of the AMDB template use 
  # the following format to output objects:
  #   report<<"object \n"<<object<<endl;
  #
  # The part in quotations becomes the list name.
  # Created By Steven Martell
  options(warn=-1)  #Suppress the NA message in the coercion to double
  
  
  ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
  idx=sapply(as.double(ifile),is.na)
  vnam=ifile[idx] #list names
  nv=length(vnam) #number of objects
  A=list()
  ir=0
  for(i in 1:nv)
  {
    ir=match(vnam[i],ifile)
    if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
    dum=NA
    if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
    if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))

    if(is.numeric(dum))#Logical test to ensure dealing with numbers
    {
      A[[vnam[i]]]=dum
    }
  }
  options(warn=0)
  
  return(A)
}

read.psv <- function(fn, nsamples=10000) {
  #This function reads the binary output from ADMB
  #-mcsave command line option.
  #fn = paste(ifile,'.psv',sep='')
  filen <- file(fn, "rb")
  nopar <- readBin(filen, what = integer(), n = 1)
  mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
  mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
  close(filen)
  return(mcmc)
}

# A simple function for creating transparent colors
# Author: Nathan Stephens (hacks package)
colr <- function(col.pal=1,a=1) {
    col.rgb<-col2rgb(col.pal)/255
    rgb(t(col.rgb),alpha=a)
}