### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(fastGHQuad)
library(flatlandr)
library(gmm)
library(Rearrangement)
library(mgcv)
library(neuralnet)
library(AER)
library(fixest)
library(lhs)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) %>%
  filter(designation !="none", !is.na(age_hh), !is.na(highest_ed), !is.na(inc_wave1),
         !is.na(imputed_bufferstock), food_consumption>0, !is.na(total_income), total_income>0) %>%
  mutate(M_avg = NA,
         total_income_hat = NA,
         vdc_id = as.factor(paste(district, vdc, ward, sep = "_")))

### purge hh heterogeneity and life cycle
avghhdat = data.frame("age_hh" = mean(df.hh$age_hh),
                      "hhmembers" = mean(df.hh$hhmembers),
                      "under18" = mean(df.hh$under18),
                      "over65" = mean(df.hh$over65),
                      "highest_ed" = mean(df.hh$highest_ed))

### income a bit tricky since don't observe previous income in wave 1
increg <- lm(log(total_income) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
               poly(highest_ed, 2), data = df.hh, weights = wt_hh)
df.hh$incresid <- increg$residuals
df.hh$total_income_hat <- exp(predict(increg, avghhdat) + increg$residuals)

### Income reg based on unaffected districts
#df.hhnone <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) %>%
#  filter(designation =="none", !is.na(age_hh), !is.na(highest_ed), !is.na(inc_wave1),
#         !is.na(imputed_bufferstock), food_consumption>0, !is.na(total_income), total_income>0) %>%
#  mutate(M_avg = NA,
#         total_income_hat = NA,
#         vdc_id = as.factor(paste(district, vdc, ward, sep = "_")))
#
#incregnone <- lm(log(total_income) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
#               poly(highest_ed, 2), data = df.hhnone, weights = wt_hh)
#df.hhnone$incresid <- incregnone$residuals
#df.hhnone$total_income_hat <- exp(predict(incregnone, avghhdat) + incregnone$residuals)
#stargazer(increg, incregnone, type = "text")

## Expected income regression - can include variables endogenous to quake bc estimating expected Y post quake (pre-aid)
einc <- feols(incresid ~ femalehh + poly(inc_wave1, 3) + poly(landvalue, 2) | vdc_id + caste_recode,
              data = df.hh, weights = df.hh$wt_hh, subset = df.hh$wave!=1)
df.hh$M_avg[df.hh$wave!=1] <- exp(predict(einc, newdata = data.frame("landvalue" = df.hh$landvalue[df.hh$wave!=1],
                                                                     "inc_wave1" = df.hh$inc_wave1[df.hh$wave!=1],
                                                                     "femalehh" = df.hh$femalehh[df.hh$wave!=1],
                                                                     "vdc_id" = df.hh$vdc_id[df.hh$wave!=1],
                                                                     "caste_recode" = df.hh$caste_recode[df.hh$wave!=1])) + predict(increg, avghhdat))

### use average estimate for waves 2 and 3 for all 3 waves - should be almost same anyway
df.avginc <- df.hh %>%
  group_by(hhid) %>% 
  summarize(M_avg = mean(M_avg, na.rm = TRUE)) 
df.hh$M_avg = NULL
df.hh <- merge(df.hh, df.avginc, by = "hhid")
#ggplot(df.hh) + geom_point(aes(x = log(M_avg), y = log(total_income)))
#ggplot(df.hh) + geom_point(aes(x = log(total_income_hat), y = log(M_avg)))

### other variables more straightforward

foodreg <- lm(log(food_consumption) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2), data = df.hh, weights = wt_hh)
df.hh$food_consumption_hat <- exp(predict(foodreg, newdata = avghhdat) + foodreg$residuals)

homereg <- tobit(home_value ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                   poly(highest_ed, 2), data = df.hh, weights = wt_hh)
df.hh$home_value_hat <- pmax(0, predict(homereg, newdata = avghhdat) + resid(homereg))

buffreg <- lm(ihs(imputed_bufferstock) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2), data = df.hh, weights = wt_hh)
buffhat <- predict(buffreg, newdata = avghhdat)
df.hh$imputed_bufferstock_hat <- .5*(exp(buffhat + buffreg$residuals) - exp(-1*(buffhat + buffreg$residuals)))

investreg <- tobit(home_investment ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                     poly(highest_ed, 2), data = df.hh, weights = wt_hh)
df.hh$home_investment_hat <- pmax(0, predict(investreg, avghhdat) + resid(investreg))

df.adj <- df.hh %>%
  select(hhid, wave, wt_hh, food_consumption_hat, home_value_hat, home_investment_hat, 
         imputed_bufferstock_hat, total_income_hat, quake_aid, M_avg, max_interest, cut_food, dissaving,
         dist_2_seg13, quake_losses, received_aid, wt_hh, gorkha_hh, NGO_transfers, designation) %>%
  mutate(food_consumption = food_consumption_hat/M_avg,
         home_value = home_value_hat/M_avg,
         home_investment = home_investment_hat/M_avg,
         imputed_bufferstock = imputed_bufferstock_hat/M_avg,
         quake_aid = quake_aid/M_avg,
         total_income = total_income_hat/M_avg) %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave),
         lag_x = lag(imputed_bufferstock, order_by = wave),
         lag_a = lag(quake_aid, order_by = wave),
         lag_y = lag(total_income, order_by = wave),
         lag_c = lag(food_consumption, order_by = wave),
         lag_i = lag(home_investment, order_by = wave)) 
saveRDS(df.adj, paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))

### generate grid points using latin hypercube sampling prior to running all jobs
set.seed(99); iterations <- as.data.frame(improvedLHS(5000, 9))
colnames(iterations) <- c("gamma", "beta", "R", "cbar", "hbar", 
                          "lambda", "sigma", "alpha", "delta")
mins <- c(1.01, .6, 1.01, .01, .01, .01, .01, .1, .6)
ranges <- c(8.99, .39, .4, .89, .89, 5.99, .99, .8, .39)
iterations <- iterations*rep(ranges, each = nrow(iterations)) + rep(mins, each = nrow(iterations))
iterations = filter(iterations, R*beta<1, cbar < alpha, hbar < 1-alpha)
saveRDS(iterations, paste(getwd(), "/data/model_output/iterations_raw.rds", sep = ""))
