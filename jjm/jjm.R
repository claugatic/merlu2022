#
library(ggplot2)
library(dplyr)
library(adnuts)
#install.packages("devtools")
#devtools::install_github("SPRFMO/jjmR")
library(jjmR)

#m1 <- runit("m1",output="results")

m1 <- runit("m1",pdf=TRUE,portrait=F,est=TRUE,exec="./jjms")

tidy_jjm <- tidy_JJM(m1)

index_fits <- tidy_jjm$index_fits

index_fits %>% 
  ggplot() + 
  geom_pointrange(aes(year, observed_ind, ymin = observed_ind - 1.96 * observed_se, ymax =  observed_ind + 1.96 * observed_se), alpha = 0.5) +
  geom_path(aes(year, pred_ind, color = model)) + 
  facet_wrap(~ fleet_name, scales = "free_y") + 
  scale_x_continuous(name = "Year", guide = guide_axis(n.dodge = 2)) + 
  scale_y_continuous(name = "Index Values")
  index_fits %>% 
  mutate(residual = pred_ind - observed_ind ) %>% 
  group_by(fleet_name, model) %>% 
  mutate(standardized_residual = residual / sd(residual)) %>% 
  filter(!is.na(standardized_residual)) %>% 
  ggplot() + 
  geom_hline(yintercept = 0,linetype = 2) +
  geom_col(aes(x = year, y =standardized_residual, fill =model), position = position_dodge(width = 0.5)) +
  facet_wrap(~ fleet_name, scales = "free_x") + 
  scale_x_continuous(name = "Year", guide = guide_axis(n.dodge = 2)) + 
  scale_y_continuous(name = "Standardized Residuals")

plot_selectivities(get_selectivities(m1))
kobe(m1, engine = "ggplot")

age_fits <- get_age_fits(m1)
age_fits
age_fits %>% 
  filter(model == "Model 1", stock == "", year > 2000) %>% 
  pivot_longer(predicted:observed) %>% 
  ggplot() + 
  geom_density(aes(age, value, fill = name),stat = "identity", alpha = 0.5) + 
  facet_grid(year~fleet_name)


recruits <- get_recruits(m1)


recruits %>% 
  ggplot() + 
  geom_ribbon(aes(year, ymin = lower_recruits, ymax = upper_recruits, fill = stock),alpha = 0.5) + 
  geom_line(aes(year, recruits, color = stock)) + 
  facet_wrap(~model)

fishing_mortality <- get_fishing_mortality(m1)

fishing_mortality %>% 
  ggplot(aes(year, mortality, color = age, group = age)) + 
  geom_line() + 
  facet_grid(model~stock, scales = "free_y") 






m <- './jjms'
# Directory  
fishing_mortality <- get_fishing_mortality(m1)

fishing_mortality %>% 
  ggplot(aes(year, mortality, color = age, group = age)) + 
  geom_line() + 
  facet_grid(model~stock, scales = "free_y") 


#---------------ADNUTS------------------------- 
m <- 'jjms'
d <- 'jjm'
setwd(d)
system(paste(m, '-nox -iprint 200 -hbf 1'))
setwd('..')
iter <- 2000 # maybe too many...depends are number cores...I used 8...
chains=8
fit.mle <- sample_nuts(model=m, path=d, iter=iter, warmup=iter/4, 
                   chains=chains, cores=chains, control=list(adapt_delta=.9,metric='mle'))
saveRDS(fit.mle,"jjm/fit.mle.RDS")
#fit.mle <- readRDS(here("fit.mle.RDS"))
adnuts::pairs_admb(fit.mle, pars=1:6, order='slow')
launch_shinyadmb(fit.mle)
# now do mceval
setwd(d)
system(paste(m, '-mceval'))
setwd('..')
# Read in some mcmc results
library(tidyverse)
library(ggthemes)
mcldf <- read_csv("mclike.csv")
mcldf
mnssb <- mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% group_by(type) %>% summarise(mean_SSB=mean(SSB)) 
mnssb
bin <- seq(5000,12000,by=500); bin
mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% mutate(bin=cut(SSB, breaks=bin, labels=FALSE)) %>% group_by(type) %>% mutate(NLL=value-mean(value)) %>%
ungroup() %>% group_by(type,bin) %>% summarise(NLL=mean(NLL))%>%
ggplot(aes(x=bin,y=NLL,color=type)) + theme_few()+ geom_line() #+ facet_wrap(type~.) 

mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% group_by(type) %>% mutate(NLL=value-mean(value)) %>%
ggplot(aes(x=SSB,y=NLL,color=type)) + theme_few()+ geom_point(alpha=.2,size=.2) + facet_wrap(type~.,scales="free") + stat_smooth()

mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% group_by(type) %>% mutate(NLL=value-mean(value)) %>%
ggplot(aes(x=SSB,y=NLL,color=type)) + theme_few()+ geom_point(alpha=.2) + facet_wrap(type~.,scales="free") + stat_smooth()

#mcldf %>% filter(stock==1,type!="ind_len",type!="priors") %>% mutate(bin=cut(SSB, breaks=bin, labels=FALSE)) %>% group_by(type) %>% mutate(NLL=value-min(value)) %>%
tmp <- mcldf %>% filter(stock==1,type!="ind_len",type!="priors") #%>% mutate(SSB=arules::discretize(SSB, breaks = 3, labels = c("Low","Medium","High"))
tmp$SSB <- arules::discretize(tmp$SSB, breaks = 5, labels = c("Low","Med-low","Medium","Med-high","High"))
tmp %>% group_by(type) %>% mutate(NLL=value-min(value)) %>% group_by(type,SSB) %>% summarise(NLL=mean(NLL)) %>%
ggplot(aes(x=SSB,y=NLL,color=type)) + theme_few()+ geom_point(size=2) + facet_wrap(type~.,scales="free") 
# mutate(bin=cut(SSB, breaks=bin, labels=FALSE)) 

mcdf <- read.table("mceval.rep",header=TRUE)
mcdf <- as.tibble(mcdf)
glimpse(mcdf)
unique(mcdf$type)
unique(mcdf$Age)
mcdf %>% filter(Year>1990,type=="SSB", Age=="all_stock_1") %>% mutate(Year=as.factor(Year) ) %>%
  ggplot(aes(x=Year,y=value)) + geom_violin(color="salmon",fill="salmon") + ggthemes::theme_few() 

mcdf %>% filter(Year>1990,type=="Recruits", Age=="Age_1_stock_1") %>% mutate(Year=as.factor(Year) ) %>%
  ggplot(aes(x=Year,y=value)) + geom_violin(color="salmon",fill="salmon") + ggthemes::theme_few() 

mcdf %>% filter(Year>1990,type=="Depletion", Age=="all_stock_1") %>% mutate(Year=as.factor(Year) ) %>%
  ggplot(aes(x=value)) + geom_density(color="salmon",fill="salmon") + ggthemes::theme_few() 
  



pairs_admb(fit) # modified pairs just for ADMB fits like this
## Can also use ShinyStan (make sure to exit it)
plot_sampler_params(fit.mle)                # NUTS adaptation
## Compare MLE and posterior of marginals. See help for
## recommendation for creating multipage PDF for high dimensional
## parameters.
plot_marginals(fit)

### ------------------------------------------------------------
### Extracting posterior samples
## Get post-warmup, merged chains into data frame (these two are
## identical)
str(as.data.frame(fit))
str(extract_samples(fit))
## If you want it in list form, e.g., to put into coda package
str(extract_samples(fit, as.list=TRUE))
## If you want to see warmup samples, and the log-posterior (lp__ column)
str(extract_samples(fit, inc_warmup=TRUE, inc_lp=TRUE))

## Remove folder
unlink(path, TRUE)