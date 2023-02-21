#                                                                                                  #
#                                                                                                  #          
#                                                                                                  #
#   Does GLP-1 receptor agonist liraglutide alter food-related sensory pleasure 
#              in patients with obesity? A randomized controlled trial                             #
#                                                                                                  #
#                                                                                                  #
#                   David Munoz Tord                                                               #
#                   Eva R Pool                                                                     #
# 

#                                        
# created  by D.M.T                                                                                #
# modified by E.R.P based on B.M advice



#-------------------------------------------------------------------------------------------------------------------------------------------
#                                       PRELIMINARY STUFF
#-------------------------------------------------------------------------------------------------------------------------------------------

# clear variables in environment
rm(list = ls())

# load libraries
if(!require(pacman)) {
  install.packages('pacman')
  library(pacman)
}

pacman::p_load(lme4, lmerTest, rstudioapi, ggplot2, grid, gridExtra, tidyverse, plyr,
               multcomp, afex, MBESS, WebPower, emmeans, Rmisc, JWileymisc, 
               kableExtra,BayesFactor,cmdstanr, brms,devtools,bayestestR,logspline)


pacman::p_load(tinylabels, apaTables, MBESS, afex, ggplot2, 
               ggpubr, Rmisc, emmeans, tidyr, BayesFactor, 
               bayestestR, devtools, lspline, kableExtra, sjPlot, 
               knitr, XML, rlist, janitor, optimx, ggthemes, dplyr, 
               JWileymisc, corrplot, caret, intmed, stringr, cowplot, 
               pander, psych, cmdstanr, brms, tidybayes, tibble, pbkrtest, scales)

# get tool
devtools::source_gist("2a1bb0133ff568cbe28d", 
                      filename = "geom_flat_violin.R")

devtools::source_gist("383aa93ffa161665c0dca1103ef73d9d", 
                      filename = "effect_CI.R")


# set path
current_dir <- dirname(getActiveDocumentContext()$path)
setwd(current_dir)
setwd('..')
home_path <- getwd()
setwd(home_path)


set_cmdstan_path(path = '/opt/anaconda3/envs/stan/bin/cmdstan')

# set working directory for outputs
analysis_path <- dirname(getActiveDocumentContext()$path)
figures_path  <- file.path(analysis_path, 'figures')
output_path <-  file.path(analysis_path, "outputs")

# source utilities
source(file.path(analysis_path,'R','utils.R'), echo=F)# useful functions

# open databases
BEHAV  <- read.delim(file.path(home_path, 'behavioral','databases',"OBIWAN_HEDONIC.txt"), header = T, sep ='') # read in dataset

#BEHAV <- read.csv(file.path(home_path,'behavioral','databases','HEDONIC.csv'), header = T, sep =',') # read in dataset
METAB <- read.csv(file.path(home_path,'behavioral','databases','METABOLIC.csv'), header = T, sep =',') # read in dataset

NEEDS <- read.csv(file.path(home_path,'behavioral','databases','internal.csv'), header = T, sep =',') # read in dataset


METAB2 <- read.csv(file.path(home_path,'behavioral','databases','Full_metabolic.csv'), header = T, sep =',') # read in dataset


# set theme for plots
averaged_theme <- theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 12, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.2, color="lightgrey"),
        panel.grid.minor = element_blank(),
        legend.text  = element_text(size =  14),
        legend.title = element_text(size =  14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =  20),
        axis.line = element_line(size = 0.5))
#panel.border = element_blank())



timeline_theme <- theme_bw(base_size = 20, base_family = "Helvetica")+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        strip.background = element_rect(color="white", fill="white", linetype="solid"),
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(size=.2, color="lightgrey"),
        panel.grid.minor = element_blank(),
        legend.text  = element_text(size =  14),
        legend.title = element_text(size =  14),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size =  20),
        axis.line = element_line(size = 0.5))
#panel.border = element_blank())

pal = viridis::inferno(n=5)


# analysis setting

my_control = lmerControl(optimizer ='optimx', optCtrl=list(method='nlminb'))

options(contrasts = c("contr.sum", "contr.poly")) # this is very important

#-------------------------------------------------------------------------------------------------------------------------------------------
#                                                        PREPROCESSING
#-------------------------------------------------------------------------------------------------------------------------------------------

METAB2 = METAB2[order(METAB2$id, decreasing = TRUE), ] 
METAB = METAB[order(METAB$id, decreasing = TRUE), ] 


# foor metabolical
METAB2$intervention <- METAB$intervention
METAB2$intervention <- dplyr::recode(METAB2$intervention, "0" = "Placebo", "1" = "Liraglutide" )



# subset only obsese
BEHAV = subset (BEHAV, group != "control")


#keep uncentered for descriptive stats
METAB$ageF = METAB$age;  

# reverse BMI_diff so it in terms of actual weight loss and not weight gain
METAB$weightLoss = -(METAB$BMI_diff); 
METAB$bmi1 = METAB$BMI_V1; METAB$dif_interv = METAB$Date_diff 

#center covariates 
biomed <- c('age','Date_diff','BW_diff','BMI_diff','WC_diff','Insulin_diff','X2.AG_diff','reelin_diff','MCP_diff','TNFalpha_diff','GLP_diff','Ghrelin_diff','Glu_diff','HOMA_IR_diff', 'AEA_diff', 'OEA_diff', 'PEA_diff', 'BMI_V1')
METAB = METAB %>% group_by %>% mutate_at(biomed, scale)

#remove outliers from biomedical (+- 3 SD)
df_dict <- data.frame(variable = biomed, out_low = rep(-3,length(biomed)),  out_high = rep(3,length(biomed)))
for (var in df_dict$variable) {
  METAB[[var]][METAB[[var]] < df_dict[df_dict$variable == var, ]$out_low | METAB[[var]] > df_dict[df_dict$variable == var, ]$out_high] <- NaN}


ALL = join(BEHAV, METAB, by = "id")


# exclude participant that did not do two sessions
ALL = subset (ALL, id != "201" &
                id != "208" &
                id != "210" &
                id != "214" &
                id != "214" &
                id != "216" &
                id != "219" &
                id != "222" &
                id != "223" &
                id != "233" &
                id != "234" &
                id != "240" &
                id != "245" &
                id != "247" &
                id != "249" &
                id != "258" &
                id != "263" &
                id != "267"
)


# remove participant that do not have fmri data
ALL = subset (ALL, id != "226" &
                id != "228" &
                id != "234" &
                id != "242"
)



fac <- c("id", "trial", "condition", "session", "intervention", "gender"); ALL[fac] <- lapply(ALL[fac], factor)


# recode to be interpretable
ALL$gender  <- dplyr::recode(ALL$gender, "0" = "Female", "1" = "Male" )
ALL$intervention <- dplyr::recode(ALL$intervention, "0" = "Placebo", "1" = "Liraglutide" )


# create averaged databases

HED.means <- aggregate(ALL$perceived_liking, by = list(ALL$id, ALL$condition, ALL$session, ALL$intervention), FUN='mean', na.rm = T) # extract means
colnames(HED.means) <- c('id','condition','session', 'intervention','perceived_liking')

HED.trial <- aggregate(ALL$perceived_liking, by = list(ALL$id, ALL$trialxcondition, ALL$condition, ALL$session, ALL$intervention), FUN='mean') # extract means
colnames(HED.trial) <- c('id',"trialxcondition",'condition',"session","intervention","perceived_liking")


#-------------------------------------------------------------------------------------------------------------------------------------------
#                                                        DEMOGRAFIC DATA
#-------------------------------------------------------------------------------------------------------------------------------------------

df.demo <- aggregate(ALL$bmi1, by = list(ALL$id, ALL$ageF, ALL$BW_V1, ALL$WC_V1, ALL$gender, ALL$intervention), FUN='mean', na.rm = T) # extract means
colnames(df.demo) <- c('id','age',"body weight", "waist circumference",'gender', 'intervention','BMI')


table_demo <- egltable(c("BMI", "age", "body weight", "waist circumference", "gender"), 
                       g = "intervention", data = df.demo, strict = FALSE) %>%
  kbl(caption ="Demografics of Trial Population Before Intervetion", digits = 2) %>%
  kable_styling(latex_options = "HOLD_position", position = "center", full_width = F) %>%
  row_spec(0,bold=T,align='c')


pdf(file.path(figures_path,'Table_Demografics.pdf'))

dev.off()



#-------------------------------------------------------------------------------------------------------------------------------------------
#                                                        METABOLIC DATA
#-------------------------------------------------------------------------------------------------------------------------------------------




# exclude participant that did not do two sessions
METAB2  = subset (METAB2 , id != "201" & id != "200" &
                id != "208" &
                id != "210" &
                id != "214" &
                id != "214" &
                id != "216" &
                id != "219" &
                id != "222" &
                id != "223" &
                id != "233" &
                id != "234" &
                id != "240" &
                id != "245" &
                id != "247" &
                id != "249" &
                id != "258" &
                id != "263" &
                id != "267")


# remove participant that do not have fmri data
METAB2 = subset (METAB2 , id != "226" &  id != "228" &  id != "234" & id != "242" & id != "212" & id != "264")




meta.subscale <- c("id","BMI_diff", "BW_diff", "WC_diff", "intervention", "Insulin_diff", "Glu_diff",
                   "SPB_diff","DBP_diff","HR_diff", "ChTot_diff","HDL_diff", "LDL_diff","TG_diff",
                   "HBA1c_diff")
df.meta <- METAB2[meta.subscale]


colnames(df.meta) <- c('id','BMI',"body weight", "waist circumference", 'intervention',
                        'insuline','gluca','systolic blood pressure', 'diastolic blood pressure',
                        'heart rate', 'total cholesterol',# missing leptin and adiopecin # restini #obestatin
                        'HDL cholesterol', 'LDL cholesterol','Triglycerides','HbA1c')

meta.variables <- c('BMI',"body weight", "waist circumference",
                    'insuline','gluca','systolic blood pressure', 'diastolic blood pressure',
                    'heart rate', 'total cholesterol',# missing leptin and adiopecin # restini #obestatin
                    'HDL cholesterol', 'LDL cholesterol','Triglycerides','HbA1c')


table_medic <- egltable(meta.variables, g = "intervention", data = df.meta, strict = FALSE) %>%
  kbl(caption ="Change in secondary end points from baseline to 16-week follow-up", digits = 2) %>%
  kable_styling(latex_options = "HOLD_position", position = "center", full_width = F) %>%
  row_spec(0,bold=T,align='c')


pdf(file.path(figures_path,'Table_Medical.pdf'))

dev.off()





#-------------------------------------------------------------------------------------------------------------------------------------------
#                                                        WEIGHT LOSS
#-------------------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------  define database
df.weight <- aggregate(ALL$BMI_diff, by = list(ALL$id, ALL$gender, ALL$age, ALL$intervention), FUN='mean', na.rm = T) # extract means
colnames(df.weight) <- c('id','gender','age', 'intervention','BMI_diff')
df.weight$gender <- ifelse(df.weight$gender=="Male",-1,1)  #SUM-CODING HERE!
df.weight$intervention <- factor(ifelse(df.weight$intervention=="Placebo",-1,1))  #SUM-CODING HERE!
df.weight$BMI_diff_z <- scale (df.weight$BMI_diff)  #SUM-CODING HERE!


#---------------------------------------------------- define statstical Model
mf = formula(BMI_diff_z ~ intervention + gender + age  + (1|id))

#---------------------------------------------------- run frequenstistic stat
fmod_weight = lm(update(mf, ~.- (1|id)) , data = df.weight)

anova(fmod_weight)

summary(fmod_weight)

confint(fmod_weight, level = 0.95, method = "Wald") 

# ----- visualize assumptions check
plot(fitted(fmod_weight),residuals(fmod_weight)) 
qqnorm(residuals(fmod_weight))
hist(residuals(fmod_weight))

#---------------------------------------------------- run baysian stat
niter = 5000; warm = 1000; chains = 4; cores = 4; nsim = 40000 

bmod_weight = brm(BMI_diff_z ~ intervention + gender + age, data=df.weight, family = gaussian, 
                  prior = c(prior(normal(0,1), class = "b")), 
                  sample_prior=TRUE, chains=4,
                  iter=niter, warmup=1000, seed=123, backend="cmdstanr",
                  control = list(adapt_delta = 0.99))
bmod_weight



# 1) Generic informative prior around 0 for fixed effects 
# 2) we need to sample priors and save parameters for computing BF 
# 3) larger step size 

# save model so that we do not have to run it evertime
save (bmod_weight, file = file.path(analysis_path, 'Bmod_weight.R'))



# extract BF
full_tab_weight = describe_posterior(bmod_weight, estimate = "median", 
                                     dispersion = T, ci = .9, ci_method = "hdi", 
                                     bf_prior = bmod_weight, diagnostic = "Rhat",  
                                     test = c("p_direction", "bf"))

# bayesian interpretation of the results (for intervetion)
report = capture.output(sexit(bmod_weight, ci=.9))
report[4]


# visulize posterior and chain sampling
param = mcmc_plot(object = bmod_weight, pars =c("b_.*en", "b_.*ge"), type ="areas")
trace = mcmc_plot(object = bmod_weight, pars =c("b_.*en", "b_.*ge"), type ="trace")



# ----- visualize assumptions check
var_group = pp_check(bmod_weight, type = "stat_grouped", group = "intervention", binwidth = 0.1, nsamples = NULL) #equality of variance between groups
rep_fit   = pp_check(bmod_weight, nsamples = 100) # check response fit
error     = pp_check(bmod_weight, type ="error_scatter_avg", nsamples = NULL) # check good alignment between model and data, and no obvious pattern to the types of errors we are getting.

diagWEIGHT <- ggarrange(param, var_group, rep_fit, error, ncol = 2, nrow = 2)

# save plot
pdf(file.path(figures_path,'BLMX_Weight_Loss_checks.pdf'))
print(diagWEIGHT)
dev.off()

# ------------------------------------------------- PLOT  ---------------------------------------------------------------------------------

# non normalized body weight for plots
ALL$BMI_Diff_raw = ALL$bmi1 - ALL$BMI_V10 
df.weight <- aggregate(ALL$BMI_Diff_raw, by = list(ALL$id, ALL$gender, ALL$age, ALL$intervention), FUN='mean', na.rm = T) # extract means
colnames(df.weight) <- c('id','gender','age', 'intervention','BMI_Diff_raw')

# compute mean and se
weight.pp <- summarySE(df.weight,
                       measurevar = c("BMI_Diff_raw"),
                       groupvars = "intervention")

# plot
p_weight = ggplot(df.weight, aes(x = intervention, y = BMI_Diff_raw, fill = intervention, color = intervention)) +
  geom_abline(slope=0, intercept=0, linetype = "dashed", color = "black") +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .1, aes(fill = intervention, color = NA), color = NA) +
  geom_point(aes(group = id, y = BMI_Diff_raw), alpha = .3, position = position_dodge(0.2), size = 0.5) +
  geom_crossbar(data = weight.pp, aes(y = BMI_Diff_raw, ymin=BMI_Diff_raw-se, ymax=BMI_Diff_raw+se), width = 0.5 , alpha = 0.1) +
  ylab('Weight Loss [BMI post - BMI pre]') +
  xlab('Intervention') +
  scale_fill_manual(values=c("Placebo"= pal[1], "Liraglutide"=pal[3]), guide = 'none') +
  scale_color_manual(values=c("Placebo"=pal[1], "Liraglutide"=pal[3]), guide = 'none') +
  theme_bw()

# save plot
p_weight + timeline_theme

pdf(file.path(figures_path,'Figure_Weight_Loss_raw.pdf'))
print(p_weight )
dev.off()


#-------------------------------------------------------------------------------------------------------------------------------------------
#                                                        TEST BEHAVIORAL DATA
#-------------------------------------------------------------------------------------------------------------------------------------------



# ---------------------------------------------------  define database

HED.aov     <- aov_car(perceived_liking ~ trialxcondition + Error (id/trialxcondition), data = HED.trial, anova_table = list(correction = "GG", es = "pes"))
HED.emm    <- emmeans(HED.aov, ~ trialxcondition , model = "multivariate")
HED.linear  <- contrast(HED.emm, "poly") # look only at the linear constrast here
coef(HED.linear)

HED.trial$saturation = scales::rescale(as.numeric(HED.trial$trialxcondition), to = c(19,-19))
HED.trial$perceived_liking_z = scale(HED.trial$perceived_liking) # so that the prior are scaled correctly
HED.trial$intervention_sc <- factor(ifelse(HED.trial$intervention=="Placebo",-1,1))  #SUM-CODING HERE!
HED.trial$condition_sc <- factor(ifelse(HED.trial$condition=="Empty",-1,1))  #SUM-CODING HERE!
HED.trial$session_sc <- factor(ifelse(HED.trial$session=="second",-1,1))  #SUM-CODING HERE!


#---------------------------------------------------- define statstical Model
mf = formula(perceived_liking_z ~ condition_sc*session_sc*intervention_sc*saturation + (condition_sc*session_sc*saturation|id))



#---------------------------------------------------- run frequenstistic stat
fmod_like = lmer(mf, data = HED.trial, control = my_control)
anova(fmod_like,type=2)
summary(fmod_like)
confint(fmod_like, level = 0.95, method = "Wald") 

# ----- visualize assumptions check
plot(fitted(fmod_like),residuals(fmod_like)) 
qqnorm(residuals(fmod_like))
hist(residuals(fmod_like))
simulateResiduals(fittedModel=fmod_like, n=1000, plot=TRUE)
plot(density(HED.trial$perceived_liking))


#---------------------------------------------------- run baysian stat
niter = 5000; warm = 1000; chains = 4; cores = 4; nsim = 40000 

my_stanvars <- stanvar(scode = stan_funs, block = "functions")

bmod_like = brm(bf(perceived_liking_z ~ condition_sc*session_sc*intervention_sc*scale(saturation) + (condition_sc*session_sc*scale(saturation)|id), hu~1),
                family = hurdle_gaussian, 
                stanvars = stanvars, 
                data=HED.trial, 
                prior =  c(prior(normal(0,1), class="b", coef=""),prior(student_t(3,0,5), class="sd")),
                sample_prior=TRUE, 
                chains = 4,  
                iter = 2000, 
                warmup = 500, 
                seed = 123, 
                backend = "cmdstan", 
                control = list(adapt_delta = 0.99)) 


#1) Custom gaussian hurdle because of value inflation at scale = 50 (neutral value selected by default on control solution)
#2) Generic weakly informative prior around 0 for fixed effects 
#3) we need to sample priors and save parameters for computing BF #this one is a big longer, so seat tight 
#4) increased step size

save (bmod_like, file = file.path(analysis_path, 'Bmod_like.R'))

full_tab_like = describe_posterior(bmod_like, estimate = "median", dispersion = T,
                                   ci = .9, ci_method = "hdi",
                                   bf_prior = bmod_like, diagnostic = "Rhat",
                                   test = c("p_direction", "bf"))


# bayesian interpretation of the results (for intervetion)
report = capture.output(sexit(bmod_like, ci=.9))
report[4]


# visulize posterior and chain sampling
param = mcmc_plot(object = bmod_like, pars =c("b_.*en", "b_.*ge"), type ="areas")
trace = mcmc_plot(object = bmod_like, pars =c("b_.*en", "b_.*ge"), type ="trace")



# ----- visualize assumptions check
var_group = pp_check(bmod_like, type = "stat_grouped", group = "intervention", binwidth = 0.1, nsamples = NULL) #equality of variance between groups
rep_fit   = pp_check(bmod_like, nsamples = 100) # check response fit
error     = pp_check(bmod_like, type ="error_scatter_avg", nsamples = NULL) # check good alignment between model and data, and no obvious pattern to the types of errors we are getting.

diagLIKE <- ggarrange(param, var_group, rep_fit, error, ncol = 2, nrow = 2)

# save plot
pdf(file.path(figures_path,'BLMX_LIKE_checks.pdf'))
print(diagLIKE)
dev.off()


# ------------------------------------------------- PLOT 1 ---------------------------------------------------------------------------------
HED.p <- summarySEwithin(HED.trial,
                         measurevar = c("perceived_liking"),
                         withinvars = c("trialxcondition","condition","session"),
                         betweenvars = "intervention",
                         idvar = "id")

# rename variables
HED.p$condition  <- dplyr::recode(HED.p$condition, "Empty" = "Neutral", "MilkShake" = "MilkShake" )
HED.p$session <- dplyr::recode(HED.p$session, "third" = "Post", "second" = "Pre" )



p_behav_time <- ggplot(HED.p , aes(x = as.numeric(trialxcondition), y = perceived_liking,
                                   color = condition, 
                                   fill  = condition))+
  facet_grid(intervention~session) +
  geom_abline(slope=0, intercept=50, linetype = "dashed", color = "black") +
  geom_line(alpha = .7, size = 1, show.legend = T) +
  geom_ribbon(aes(ymax = perceived_liking + se, ymin = perceived_liking - se, fill = condition, color =NA),color = NA,  alpha=0.4) + 
  geom_point() +
  ylab('Pleasantness')+
  xlab('Trial') +
  scale_y_continuous(expand = c(0, 0), breaks = c(seq.int(0,100, by = 10)), limits = c(0,100)) +
  scale_fill_manual(values=c("Neutral"= pal[1], "MilkShake"=pal[4]), guide = 'none') +
  scale_color_manual(values=c("Neutral"=pal[1], "MilkShake"=pal[4]), guide = 'none') +
  theme_bw()


p_behav_time = p_behav_time + timeline_theme

pdf(file.path(figures_path,'Figure_Behav_time.pdf'))
print(p_behav_time )
dev.off()


# ------------------------------------------------- PLOT 2 ---------------------------------------------------------------------------------

# rename variables
HED.means$condition  <- dplyr::recode(HED.means$condition, "Empty" = "Neutral", "MilkShake" = "MilkShake" )
HED.means$session <- dplyr::recode(HED.means$session, "third" = "Post", "second" = "Pre" )


HED.pp <- summarySEwithin(HED.means,
                          measurevar = c("perceived_liking"),
                          withinvars = c("condition","session"),
                          betweenvars = "intervention",
                          idvar = "id")


p_behav_individual = ggplot(HED.means, aes(x = condition, y = perceived_liking, fill = condition, color = condition)) +
  facet_grid(intervention~factor(session, level = c("Pre","Post"))) +
  geom_abline(slope=0, intercept=50, linetype = "dashed", color = "black") +
  geom_line (aes(group = id, y = perceived_liking), alpha = .2, position = position_dodge(0.2), color = "gray") +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .1, aes(fill = condition, color = NA), color = NA) +
  geom_point(aes(group = id, y = perceived_liking), alpha = .4, position = position_dodge(0.2), size = 0.5) +
  geom_crossbar(data = HED.pp, aes(y = perceived_liking, ymin=perceived_liking-se, ymax=perceived_liking+se), width = 0.5 , alpha = 0.1) +
  ylab('Pleasantness') +
  xlab('Taste Stimulus') +
  scale_y_continuous(expand = c(0, 0), breaks = c(seq.int(0,100, by = 10)), limits = c(0,100)) +
  scale_fill_manual(values=c("Neutral"= pal[1], "MilkShake"=pal[4]), guide = 'none') +
  scale_color_manual(values=c("Neutral"=pal[1], "MilkShake"=pal[4]), guide = 'none') +
  theme_bw()

p_behav_individual = p_behav_individual + timeline_theme


pdf(file.path(figures_path,'Figure_Behav_average.pdf'))
print(p_behav_individual)
dev.off()




# -------------------------------------------------------------------------------------------------------
#                           REVISIONS
# --------------------------------------------------------------------------------------------------------


#  ------------------------- REVIEWER REQUEST 1 --------------------------------
# Run a post-hoc power analysis - this is not ideal because post-hoc powver analysis are
# basically another way to report p-value. What we can do is a sensitivity analysis explaning
# how likely it was to observe a significant effect, given your sample, and given an expected effect size 
# (see Lakens 2022, Sample Size Justification, PsyArix)


HED.sensitivity     <- aov_car(perceived_liking_z ~ condition*session*intervention*saturation + Error (id/condition*session*saturation), data = HED.trial, anova_table = list(correction = "GG", es = "pes"))
# sensitivity analysis on anova determines our studies is powered to observe f effect size of 0.22 (which can be described as a medium effect size)




# The point is an important one and we can address with more baysian modeling on the effect of interest


niter = 5000; warm = 1000; chains = 4; cores = 4; nsim = 10000 

my_stanvars <- stanvar(scode = stan_funs, block = "functions")

bmod_like_simple = brm(bf(perceived_liking_z ~ condition*session*intervention + (condition*session|id), hu~1),
                family = hurdle_gaussian, 
                stanvars = stanvars, 
                data=HED.trial, 
                prior =  c(prior(normal(0,1), class="b", coef=""),prior(student_t(3,0,5), class="sd")),
                sample_prior=TRUE, 
                chains = 4,  
                iter = 2000, 
                warmup = 500, 
                seed = 123, 
                backend = "cmdstan", 
                control = list(adapt_delta = 0.99)) 

full_tab_like_simple = describe_posterior(bmod_like_simple, estimate = "median", dispersion = T,
                                          ci = .9, ci_method = "hdi",
                                          bf_prior = bmod_like_simple, diagnostic = "Rhat",
                                          test = c("p_direction", "bf"))

param = mcmc_plot(object = bmod_like_simple, pars =c("b_condition1:session1:intervention1*"), type ="areas")




#  ------------------------- REVIEWER REQUEST 2 --------------------------------

# a plot of intervention by session 

# create averaged databases
HED.means_request1 <- aggregate(ALL$perceived_liking, by = list(ALL$id, ALL$session, ALL$intervention), FUN='mean', na.rm = T) # extract means
colnames(HED.means_request1) <- c('id','session', 'intervention','perceived_liking')
HED.means_request1$session <- dplyr::recode(HED.means_request1$session, "third" = "Post", "second" = "Pre" )


HED.pp_request1 <- summarySEwithin(HED.means_request1,
                          measurevar = c("perceived_liking"),
                          withinvars = c("session"),
                          betweenvars = "intervention",
                          idvar = "id")



p_behav_request1 = ggplot(HED.means_request1, aes(x = session, y = perceived_liking, fill = session, color = session)) +
  facet_grid(~intervention) +
  geom_abline(slope=0, intercept=50, linetype = "dashed", color = "black") +
  geom_line (aes(group = id, y = perceived_liking), alpha = .2, position = position_dodge(0.2), color = "gray") +
  geom_flat_violin(scale = "count", trim = FALSE, alpha = .1, aes(fill = session, color = NA), color = NA) +
  geom_point(aes(group = id, y = perceived_liking), alpha = .4, position = position_dodge(0.2), size = 0.5) +
  geom_crossbar(data = HED.pp_request1, aes(y = perceived_liking, ymin=perceived_liking-se, ymax=perceived_liking+se), width = 0.5 , alpha = 0.1) +
  ylab('Pleasantness') +
  xlab('Taste Stimulus') +
  scale_y_continuous(expand = c(0, 0), breaks = c(seq.int(0,100, by = 10)), limits = c(0,100)) +
  scale_fill_manual(values=c("Pre"= pal[1], "Post"=pal[4]), guide = 'none') +
  scale_color_manual(values=c("Pre"=pal[1], "Post"=pal[4]), guide = 'none') +
  theme_bw()

p_behav_request1 = p_behav_request1 + timeline_theme

# save plot
pdf(file.path(figures_path,'Plot_Revisions.pdf'))
print(p_behav_request1)
dev.off()





#  ------------------------- REVIEWER REQUEST 3 --------------------------------
# Run an analysis on the Homa-IR  rather than BMI



# ---------------------------------------------------  define database
df.homair <- aggregate(ALL$HOMA_IR_diff, by = list(ALL$id, ALL$gender, ALL$age, ALL$intervention), FUN='mean', na.rm = T) # extract means
colnames(df.homair) <- c('id','gender','age', 'intervention','HOMA_IR_diff')


#---------------------------------------------------- define statstical Model
mf = formula(HOMA_IR_diff ~ intervention + gender + age  + (1|id))

#---------------------------------------------------- run frequenstistic stat
fmod_homair = lm(update(mf, ~.- (1|id)) , data = df.homair)

anova(fmod_homair)

summary(fmod_homair)

confint(fmod_homair, level = 0.95, method = "Wald") 

# ----- visualize assumptions check
plot(fitted(fmod_homair),residuals(fmod_homair)) 
qqnorm(residuals(fmod_homair))
hist(residuals(fmod_homair))

#---------------------------------------------------- run baysian stat
niter = 5000; warm = 1000; chains = 4; cores = 4; nsim = 10000 

df.homair$gender <- ifelse(df.homair$gender==0,-1,df.homair$gender)  #SUM-CODING HERE!
df.homair$intervention <- ifelse(df.homair$intervention==0,-1,df.homair$intervention)  #SUM-CODING HERE!



bmod_homair = brm(HOMA_IR_diff ~ intervention + gender + age, data=df.homair, family = gaussian, 
                  prior = c(prior(normal(0,1), class = "b")), 
                  sample_prior=TRUE, chains=4,
                  iter=niter, warmup=1000, seed=123, backend="cmdstanr",
                  control = list(adapt_delta = 0.99))




# 1) Generic informative prior around 0 for fixed effects 
# 2) we need to sample priors and save parameters for computing BF 
# 3) larger step size 

# save model so that we do not have to run it evertime
#save (bmod_weight, file = file.path(analysis_path, 'Bmod_weight.R'))



# extract BF
full_tab_homair = describe_posterior(bmod_homair, estimate = "median", 
                                     dispersion = T, ci = .9, ci_method = "hdi", 
                                     bf_prior = bmod_homair, diagnostic = "Rhat",  
                                     test = c("p_direction", "bf"))

# bayesian interpretation of the results (for intervetion)
report = capture.output(sexit(bmod_homair, ci=.9))
report[4]


# visulize posterior and chain sampling
param = mcmc_plot(object = bmod_homair, pars =c("b_.*en", "b_.*ge"), type ="areas")
trace = mcmc_plot(object = bmod_homair, pars =c("b_.*en", "b_.*ge"), type ="trace")


# ----- visualize assumptions check
var_group = pp_check(bmod_homair, type = "stat_grouped", group = "intervention", binwidth = 0.1, nsamples = NULL) #equality of variance between groups
rep_fit   = pp_check(bmod_homair, nsamples = 100) # check response fit
error     = pp_check(bmod_homair, type ="error_scatter_avg", nsamples = NULL) # check good alignment between model and data, and no obvious pattern to the types of errors we are getting.

diagHOMA_IR <- ggarrange(param, var_group, rep_fit, error, ncol = 2, nrow = 2)

# save plot
pdf(file.path(figures_path,'BLMX_HOMA_IR_checks.pdf'))
print(diagHOMA_IR)
dev.off()



# -------------------------- REVIWER REQUEST 4 --------------------------------

# add homa ir in the analysis

HED.trial_request <- aggregate(ALL$perceived_liking, by = list(ALL$id, ALL$trialxcondition, ALL$condition, ALL$session, ALL$intervention, ALL$HOMA_IR_diff), FUN='mean') # extract means
colnames(HED.trial_request) <- c('id',"trialxcondition",'condition',"session","intervention","homa_ir","perceived_liking")

HED.aov     <- aov_car(perceived_liking ~ trialxcondition + Error (id/trialxcondition), data = HED.trial_request, anova_table = list(correction = "GG", es = "pes"))
HED.emm    <- emmeans(HED.aov, ~ trialxcondition , model = "multivariate")
HED.linear  <- contrast(HED.emm, "poly") # look only at the linear constrast here
coef(HED.linear)

HED.trial_request$saturation = scales::rescale(as.numeric(HED.trial_request$trialxcondition), to = c(19,-19))
HED.trial_request$perceived_liking_z = scale(HED.trial_request$perceived_liking) # so that the prior are scaled correctly

#---------------------------------------------------- define statstical Model

mf = formula(perceived_liking_z ~ condition*session*intervention*saturation* homa_ir+ (condition*session*saturation|id))


#---------------------------------------------------- run frequenstistic stat
fmod_like = lmer(mf, data = HED.trial_request, control = my_control)
anova(fmod_like,type=2)
summary(fmod_like)
confint(fmod_like, level = 0.95, method = "Wald") 

# ----- visualize assumptions check
plot(fitted(fmod_like),residuals(fmod_like)) 
qqnorm(residuals(fmod_like))
hist(residuals(fmod_like))
simulateResiduals(fittedModel=fmod_like, n=1000, plot=TRUE)
plot(density(HED.trial$perceived_liking))




# -------------------------- REVIEWER REQUEST 5 ----------------------------------------------------------
# correlation between change in liking and change in weight

HED.supp <- aggregate(ALL$perceived_liking, by = list(ALL$id, ALL$condition, ALL$session, ALL$intervention, ALL$BMI_diff), FUN='mean', na.rm = T) # extract means
colnames(HED.supp) <- c('id','condition','session', 'intervention','BMI_diff','perceived_liking')

HED.supp <- ddply(HED.supp, .(id), transform, liking_diff  = perceived_liking[session=="second"] - perceived_liking[session=="third"]) 

HED.supp <- subset(HED.supp, session == "third")

HED.choco <- subset(HED.supp, condition == "MilkShake")
scatterplot(HED.choco$BMI_diff, HED.choco$liking_diff)
cor.test(HED.choco$BMI_diff, HED.choco$liking_diff)

HED.neutral <- subset(HED.supp, condition == "Empty")
scatterplot(HED.neutral$BMI_diff, HED.neutral$liking_diff)
cor.test(HED.neutral$BMI_diff, HED.neutral$liking_diff)


mf = formula(perceived_liking ~ condition*session*intervention*BMI_diff+ (condition*session*intervention|id))

#---------------------------------------------------- run frequenstistic stat
fmod_like = lmer(mf, data = ALL, control = my_control)
anova(fmod_like,type=2)
summary(fmod_like)
confint(fmod_like, level = 0.95, method = "Wald") 


