library(jjmR)
library(ggplot2)
library(dplyr)
library(adnuts)

m1 <- runit("m1",output="results")

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






m <- './jjms'
# Directory  
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

pairs_admb(fit) # modified pairs just for ADMB fits like this
## Can also use ShinyStan (make sure to exit it)
## launch_shinyadmb(fit)
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