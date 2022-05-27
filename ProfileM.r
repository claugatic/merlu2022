library(r4ss)

  ## Not run: 
# note: don't run this in your main directory
# make a copy in case something goes wrong
#mydir <- "c:/CEGM/MODEL/merluza/M21/2021/4FSDLPRO/"
#c:\CEGM\MODEL\merluza\M21\E2022\m00\

mydir <- "c:/CEGM/MODEL/merluza/M21/E2022/m00/"

# the following commands related to starter.ss could be done by hand
# read starter file
starter <- SS_readstarter(file.path(mydir, 'starter.ss'))
# change control file name in the starter file
starter$ctlfile <- "control_modified.ss"
# make sure the prior likelihood is calculated
# for non-estimated quantities
starter$prior_like <- 1
# write modified starter file
SS_writestarter(starter, dir=mydir, overwrite=TRUE)

# vector of values to profile over
#h.vec <- seq(0.5,0.9,.1)
h.vec <- seq(0.1, 0.4, 0.025)
Nprofile <- length(h.vec)

# run SS_profile command
profile <- SS_profile(dir=mydir, # directory
                      # "NatM" is a subset of one of the
                      # parameter labels in control.ss_new
                      model="ss",
                      masterctlfile="control.ss_new",
                      newctlfile="control_modified.ss",
                      string="NatM_p_1_Fem_GP_1",
                      profilevec=h.vec)


# read the output files (with names like Report1.sso, Report2.sso, etc.)
profilemodels <- SSgetoutput(dirvec=mydir, keyvec=1:Nprofile)
# summarize output
profilesummary <- SSsummarize(profilemodels)

# OPTIONAL COMMANDS TO ADD MODEL WITH PROFILE PARAMETER ESTIMATED
#MLEmodel <- SS_output("C:/ss/SSv3.24l_Dec5/Simple")
#profilemodels$MLE <- MLEmodel
#profilesummary <- SSsummarize(profilemodels)
# END OPTIONAL COMMANDS

# plot profile using summary created above

#SSplotProfile(profilesummary,           # summary object
 #             profile.string = "steep", # substring of profile parameter
  #            profile.label="Stock-recruit steepness (h)") # axis label

# make timeseries plots comparing models in profile
#SSplotComparisons(profilesummary,legendlabels=paste("h =",h.vec))


## End(Not run)
#steep

dp='c:/CEGM/MODEL/merluza/M21/E2022/plotprofile/'

SSplotProfile(profilesummary, plot = TRUE, print = TRUE,
models = "all", profile.string = "NatM_p_1_Fem_GP_1",
profile.label = "M", exact = TRUE,
ylab = "Change in -log-likelihood", components = c("TOTAL", "Catch",
"Equil_catch", "Survey", "Discard", "Mean_body_wt", "Length_comp",
"Age_comp", "Size_at_age", "SizeFreq", "Morphcomp", "Tag_comp",
"Tag_negbin", "Recruitment", "InitEQ_Regime", "Forecast_Recruitment",
"Parm_priors", "Parm_softbounds", "Parm_devs", "F_Ballpark",
"Crash_Pen"), component.labels = c("Total", "Catch",
"Equilibrium catch", "Index data", "Discard", "Mean body weight",
"Length data", "Age data", "Size-at-age data", "Generalized size data",
"Morph composition data", "Tag recapture distribution",
"Tag recapture total", "Recruitment", "Initital equilibrium recruitment",
"Forecast recruitment", "Priors", "Soft bounds", "Parameter deviations",
"F Ballpark", "Crash penalty"), minfraction = 0.01,
sort.by.max.change = TRUE, col = "default", pch = "default",
lty = 1, lty.total = 1, lwd = 2, lwd.total = 3, cex = 1,
cex.total = 1.5, xlim = "default", ymax = "default", xaxs = "r",
yaxs = "r", type = "o", legend = TRUE, legendloc = "topright",
pwidth = 6.5, pheight = 5, punits = "in", res = 300,
ptsize = 10, cex.main = 1, plotdir = dp, add_cutoff = FALSE,
cutoff_prob = 0.95, verbose = TRUE)



