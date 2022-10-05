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
  group_by(hhid) %>% summarize(M_avg = mean(M_avg, na.rm = TRUE)) 
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
         imputed_bufferstock_hat, total_income_hat, quake_aid, M_avg, 
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
         lag_y = lag(total_income, order_by = wave),
         lag_c = lag(food_consumption, order_by = wave),
         lag_i = lag(home_investment, order_by = wave)) 
saveRDS(df.adj, paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))
df.adj <- df.adj[complete.cases(df.adj),]


### Grid size
xn = 50; hn = 50

### generate starting points - run 30 versions of this on cluster using dSQ 
### dsq --job-file ~/project/HPC_WBHRVS_DSQ.txt -c 20 --mem-per-cpu 2g -t 24:00:00 --mail-type ALL

#slurmseed <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
#set.seed(slurmseed); seeds = sample(c(1:1e6),9)
#set.seed(seeds[1]); gamma <- runif(500, 1.01, 10)
#set.seed(seeds[2]); beta <- runif(500, .6, .99)
#set.seed(seeds[3]); R <- runif(500, 1, 1.4)
#set.seed(seeds[4]); cbar <- runif(500, .01, .9)
#set.seed(seeds[5]); hbar <- runif(500, .01, .9)
#set.seed(seeds[6]); lambda <- runif(500, .01, 10)
#set.seed(seeds[7]); sigma <- runif(500, .01, 1)
#set.seed(seeds[8]); alpha <- runif(500, .1, .9)
#set.seed(seeds[9]); delta <- runif(500, .6, .99)

#iterations <- data.frame(gamma, beta, R, cbar, hbar,
#                         lambda, sigma, alpha, delta)
#iterations = filter(iterations, R*beta<1, cbar < alpha, hbar < 1-alpha)

gpenalty = 27946515 ## max iterations value

globalwrap <- function(t){ #chg back to ridx
  if (t[2]*t[3] > 1 | t[4] > t[8] | t[5] > 1 - t[8]){
    return(gpenalty)
  } else {
    #t = as.vector(t(iterations[ridx, ]))
    momentmat = gmmmomentmatcher(t, df.adj)
    g = colSums(momentmat); print(g)
    return(sum(g^2)) ## change this to g
  }
}

sol <- crs2lm(x0 = c(1.2, .926, 1.054, .33, .19, .08, .42, .63, .93), 
              fn = globalwrap, 
              lower = c(1.01, .6, 1, .01, .01, .01, .01, .1, .5), 
              upper = c(10, .99, 1.2, .9, .9, 10, 1.5, .99, 1),
              xtol_rel = 1e-2)
sol

#f <- lapply(c(1:nrow(iterations)), globalwrap)
#f <- do.call(rbind, f)
#colnames(f) <- paste("m", seq(1,11,1), sep = "_")
#iterations <- cbind(iterations, f)

#saveRDS(iterations, paste(getwd(), "/data/model_output/iterations", slurmseed, ".rds", sep = ""))

### after all jobs finish, run
# setwd("data/model_output")
# ifiles <- list.files(pattern = "^iter")
# ilist <- lapply(ifiles, readRDS)
# iterations <- do.call(rbind, ilist)
# saveRDS(iterations, "iterations.rds")
