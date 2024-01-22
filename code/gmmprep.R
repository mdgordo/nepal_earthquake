### packages and data
library(tidyverse)
library(nloptr)
library(fastGHQuad)
#devtools::install_github("mdgordo/flatlandr",ref="master")
library(flatlandr)
library(gmm)
library(Rearrangement)
library(fixest)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) %>%
  filter(!is.na(age_hh), !is.na(highest_ed), !is.na(religion), !is.na(home_investment), 
         !is.na(credit), !is.na(liq_savings), !is.na(tot_savings), !is.na(remittance_income), 
         !is.na(consumption), !is.na(total_income), !is.na(avg_income), !is.na(avg_remit), n_waves==3) %>%
  mutate(M_avg = NA,
         total_income_hat = NA,
         vdc_id = as.factor(paste(district, vdc, ward, sep = "_")),
         eth_id = as.factor(paste(caste_recode, religion, sep = "_")))

### purge hh heterogeneity and life cycle
avghhdat = data.frame("age_hh" = mean(df.hh$age_hh),
                      "hhmembers" = mean(df.hh$hhmembers),
                      "under18" = mean(df.hh$under18),
                      "over65" = mean(df.hh$over65),
                      "highest_ed" = mean(df.hh$highest_ed),
                      "wave" = 1)

## Expected income regression 
increg <- feols(log(total_income+1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
               poly(highest_ed, 2) + as.factor(wave) +
                 femalehh + poly(landvalue, 2) + prev_migrants + migrant_male + migrant_educated +
                 migrant_kat + migrant_india + migrant_other_asia + migrant_mideast + migrant_working +
                 migrant_earnings | eth_id + vdc_id, 
             data = df.hh, weights = df.hh$wt_hh)

df.hh$M_avg <- exp(predict(increg, newdata = data.frame("age_hh" = avghhdat$age_hh,
                                                      "hhmembers" = avghhdat$hhmembers,
                                                      "under18" = avghhdat$under18,
                                                      "over65" = avghhdat$over65,
                                                      "highest_ed" = avghhdat$highest_ed,
                                                      "wave" = avghhdat$wave,
                                                      "landvalue" = df.hh$landvalue,
                                                     "femalehh" = df.hh$femalehh,
                                                     "prev_migrants" = df.hh$prev_migrants,
                                                     "migrant_male" = df.hh$migrant_male,
                                                     "migrant_educated" = df.hh$migrant_educated,
                                                     "migrant_kat" = df.hh$migrant_kat,
                                                     "migrant_india" = df.hh$migrant_india,
                                                     "migrant_other_asia" = df.hh$migrant_other_asia,
                                                     "migrant_mideast" = df.hh$migrant_mideast,
                                                     "migrant_working" = df.hh$migrant_working,
                                                     "migrant_earnings" = df.hh$migrant_earnings,
                                                     "eth_id" = df.hh$eth_id,
                                                     "vdc_id" = df.hh$vdc_id)))
df.hh$total_income_hat <- exp(log(df.hh$M_avg) + increg$residuals)
incratio <- mean(df.hh$total_income)/mean(df.hh$total_income_hat)
df.hh$total_income_hat <- df.hh$total_income_hat*incratio

# ## checks
summary(increg)
sum(is.na(df.hh$M_avg) | df.hh$M_avg==0)
sum(is.na(df.hh$total_income_hat) | df.hh$total_income_hat==0)
## fitted values - fraction of variation excluded
ggplot() + geom_point(aes(x = log(df.hh$total_income+1), y = increg$fitted.values, color = df.hh$gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
## true vals vs adjusted vals (expected and realized)
ggplot() + geom_point(aes(x = log(1+ df.hh$total_income), y = log(1+df.hh$total_income_hat), color = df.hh$gorkha_hh))
ggplot() + geom_point(aes(x = log(1+ df.hh$total_income), y = log(df.hh$M_avg), color = df.hh$gorkha_hh))
## Expected adjusted vals vs realized adjusted vals - how much is permanent vs transitory
ggplot() + geom_point(aes(x = log(df.hh$M_avg), y = log(df.hh$total_income_hat), color = df.hh$gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

### Expected remittances + informal transfers
df.hh$remittance_income <- df.hh$remittance_income + df.hh$inf_transfers

erem <- feols(remittance_income ~ log(total_income+1) + lag_remitt + log(total_income+1):(femalehh + 
                                   prev_migrants + migrant_male + migrant_educated +  migrant_kat + migrant_india + migrant_other_asia + 
                                   migrant_mideast + migrant_working + migrant_earnings) | as.factor(hhid),
                data = df.hh, weights = df.hh$wt_hh)

### Big shock
bshock <- 2*weighted.mean(sqrt(df.hh$var_log_income), df.hh$wt_hh)

## Expected remittances
df.hh$E_remittance <- predict(erem, newdata = data.frame("total_income" = exp(df.hh$avg_log_income),
                                                         "lag_remitt" = df.hh$avg_remit,
                                                         "landvalue" = df.hh$landvalue,
                                                         "femalehh" = df.hh$femalehh,
                                                         "prev_migrants" = df.hh$prev_migrants,
                                                         "migrant_male" = df.hh$migrant_male,
                                                         "migrant_educated" = df.hh$migrant_educated,
                                                         "migrant_kat" = df.hh$migrant_kat,
                                                         "migrant_india" = df.hh$migrant_india,
                                                         "migrant_other_asia" = df.hh$migrant_other_asia,
                                                         "migrant_mideast" = df.hh$migrant_mideast,
                                                         "migrant_working" = df.hh$migrant_working,
                                                         "migrant_earnings" = df.hh$migrant_earnings,
                                                         "hhid" = df.hh$hhid))

df.hh$avail_remittance <- predict(erem, newdata = data.frame("total_income" = exp(df.hh$avg_log_income - bshock),
                                                             "lag_remitt" = df.hh$lag_remitt,
                                                             "landvalue" = df.hh$landvalue,
                                                             "femalehh" = df.hh$femalehh,
                                                             "prev_migrants" = df.hh$prev_migrants,
                                                             "migrant_male" = df.hh$migrant_male,
                                                             "migrant_educated" = df.hh$migrant_educated,
                                                             "migrant_kat" = df.hh$migrant_kat,
                                                             "migrant_india" = df.hh$migrant_india,
                                                             "migrant_other_asia" = df.hh$migrant_other_asia,
                                                             "migrant_mideast" = df.hh$migrant_mideast,
                                                             "migrant_working" = df.hh$migrant_working,
                                                             "migrant_earnings" = df.hh$migrant_earnings,
                                                             "hhid" = df.hh$hhid))

### Expected remittances - log
erem_1 <- feols(log(1+remittance_income) ~ log(total_income+1) + lag_remitt + log(total_income+1):(femalehh + 
                             prev_migrants + migrant_male + migrant_educated +  migrant_kat + migrant_india + migrant_other_asia + 
                             migrant_mideast + migrant_working + migrant_earnings) | as.factor(hhid),
               data = df.hh, weights = df.hh$wt_hh)

## Expected remittances
df.hh$E_remittance_1 <- exp(predict(erem_1, newdata = data.frame("total_income" = exp(df.hh$avg_log_income),
                                                   "lag_remitt" = df.hh$avg_remit,
                                                   "landvalue" = df.hh$landvalue,
                                                   "femalehh" = df.hh$femalehh,
                                                   "prev_migrants" = df.hh$prev_migrants,
                                                   "migrant_male" = df.hh$migrant_male,
                                                   "migrant_educated" = df.hh$migrant_educated,
                                                   "migrant_kat" = df.hh$migrant_kat,
                                                   "migrant_india" = df.hh$migrant_india,
                                                   "migrant_other_asia" = df.hh$migrant_other_asia,
                                                   "migrant_mideast" = df.hh$migrant_mideast,
                                                   "migrant_working" = df.hh$migrant_working,
                                                   "migrant_earnings" = df.hh$migrant_earnings,
                                                   "hhid" = df.hh$hhid)))

df.hh$avail_remittance_1 <- exp(predict(erem_1, newdata = data.frame("total_income" = exp(df.hh$avg_log_income - bshock),
                                                             "lag_remitt" = df.hh$lag_remitt,
                                                             "landvalue" = df.hh$landvalue,
                                                             "femalehh" = df.hh$femalehh,
                                                             "prev_migrants" = df.hh$prev_migrants,
                                                             "migrant_male" = df.hh$migrant_male,
                                                             "migrant_educated" = df.hh$migrant_educated,
                                                             "migrant_kat" = df.hh$migrant_kat,
                                                             "migrant_india" = df.hh$migrant_india,
                                                             "migrant_other_asia" = df.hh$migrant_other_asia,
                                                             "migrant_mideast" = df.hh$migrant_mideast,
                                                             "migrant_working" = df.hh$migrant_working,
                                                             "migrant_earnings" = df.hh$migrant_earnings,
                                                             "hhid" = df.hh$hhid)))

### model comparisons and checks
sum(df.hh$E_remittance<=0 | is.na(df.hh$E_remittance))
sum(df.hh$E_remittance_1<=0 | is.na(df.hh$E_remittance_1))

sum(df.hh$avail_remittance[df.hh$wave!=1]<=0 | is.na(df.hh$avail_remittance[df.hh$wave!=1]))
sum(df.hh$avail_remittance_1[df.hh$wave!=1]<=0 | is.na(df.hh$avail_remittance_1[df.hh$wave!=1]))

### How do two models compare to each other
ggplot(df.hh) +  
  geom_point(aes(x = log(E_remittance), y = log(E_remittance_1), color = log(1+remittance_income))) + 
  geom_abline(slope = 1, intercept = 0) + scale_color_viridis_c()
cor(df.hh$E_remittance, df.hh$E_remittance_1, method = "spearman")

ggplot(df.hh) +  
  geom_point(aes(x = log(avail_remittance), y = log(avail_remittance_1), color = log(1+remittance_income))) + 
  geom_abline(slope = 1, intercept = 0) + scale_color_viridis_c()
cor(df.hh$avail_remittance[df.hh$wave!=1], df.hh$avail_remittance_1[df.hh$wave!=1], method = "spearman")

etable(erem, erem_1)

### How do they compare to data
ggplot(df.hh) + 
  geom_point(aes(x = log(1+remittance_income), y = log(E_remittance), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) + 
  geom_point(aes(x = log(1+remittance_income), y = log(E_remittance_1), color = lag_remitt), alpha = .3) + 
  geom_abline(slope = 1, intercept = 0) + scale_color_viridis_c()

ggplot(df.hh) + 
  geom_point(aes(x = log(1+remittance_income), y = log(avail_remittance), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) + 
  geom_point(aes(x = log(1+remittance_income), y = log(avail_remittance_1), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

### are actual remittances greater than available remittances? 
sum(df.hh$avail_remittance < df.hh$remittance_income, na.rm = TRUE)/nrow(filter(df.hh, wave!=1))
sum(df.hh$avail_remittance_1 < df.hh$remittance_income, na.rm = TRUE)/nrow(filter(df.hh, wave!=1))
### avail less than expected? Could be ok because it's basically a debt
sum(df.hh$avail_remittance - df.hh$E_remittance <=0, na.rm = TRUE)
sum(df.hh$avail_remittance_1 - df.hh$E_remittance_1 <=0, na.rm = TRUE)

### negative obs bounded by 0 avail remittances minus expected
ggplot(df.hh) +
  geom_point(aes(x = ihs(avail_remittance_1 - E_remittance_1), y = ihs(E_remittance_1)))

### Purge expected remittances 
erempreg <- lm(log(E_remittance_1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                 poly(highest_ed, 2), 
               data = df.hh, weights = wt_hh)
df.hh$E_remittance_hat <- exp(predict(erempreg, newdata = avghhdat) + erempreg$residuals)
remratio <- mean(df.hh$E_remittance_1)/mean(df.hh$E_remittance_hat)
df.hh$E_remittance_hat <- df.hh$E_remittance_hat*remratio

summary(erempreg)
ggplot(df.hh) +
  geom_point(aes(x = log(E_remittance_1), y = erempreg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = log(E_remittance_1), y = log(E_remittance_hat), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

### Scale available remittances by ratio
df.hh$avail_remittance_hat = df.hh$avail_remittance_1 * df.hh$E_remittance_hat/df.hh$E_remittance_1

ggplot(df.hh, aes(x = log(avail_remittance_hat), y = log(avail_remittance_1))) + geom_point()

### Purge consumption, housing, default, savings variables
foodreg <- lm(log(food_consumption+1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2) + as.factor(wave), 
              data = df.hh, weights = wt_hh)
df.hh$food_consumption_hat <- exp(predict(foodreg, newdata = avghhdat) + foodreg$residuals)
foodratio <- mean(df.hh$food_consumption)/mean(df.hh$food_consumption_hat)
df.hh$food_consumption_hat <- df.hh$food_consumption_hat*foodratio

ggplot(df.hh) +
  geom_point(aes(x = log(food_consumption), y = foodreg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = log(food_consumption), y = log(food_consumption_hat), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

homereg <- lm(log(home_value+1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2) + as.factor(wave), 
              data = df.hh, weights = wt_hh)
df.hh$home_value_hat <- exp(predict(homereg, newdata = avghhdat) + homereg$residuals)
homeratio <- mean(df.hh$home_value)/mean(df.hh$home_value_hat)
df.hh$home_value_hat <- df.hh$home_value_hat*homeratio

ggplot(df.hh) +
  geom_point(aes(x = log(home_value + 1), y = homereg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = log(home_value + 1), y = log(home_value_hat), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

defreg <- lm(default ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
               poly(highest_ed, 2) + as.factor(wave), 
             data = df.hh, weights = wt_hh)
df.hh$default_hat <- predict(defreg, newdata = avghhdat) + defreg$residuals
defratio <- mean(df.hh$default)/mean(df.hh$default_hat)
df.hh$default_hat <- df.hh$default_hat*defratio

ggplot(df.hh) +
  geom_point(aes(x = default, y = defreg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = default, y = default_hat, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

savereg <- lm(log(1+tot_savings) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2) + as.factor(wave), 
              data = df.hh, weights = wt_hh)
df.hh$tot_savings_hat <- exp(predict(savereg, newdata = avghhdat) + savereg$residuals)
savratio <- mean(df.hh$tot_savings)/mean(df.hh$tot_savings_hat)
df.hh$tot_savings_hat <- df.hh$tot_savings_hat*savratio

ggplot(df.hh) +
  geom_point(aes(x = log(1+tot_savings), y = savereg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = log(1+tot_savings), y = log(tot_savings_hat), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

liqreg <- lm(log(1+liq_savings) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
               poly(highest_ed, 2) + as.factor(wave), 
             data = df.hh, weights = wt_hh)
df.hh$liq_savings_hat <- exp(predict(liqreg, newdata = avghhdat) + liqreg$residuals)
liqratio <- mean(df.hh$liq_savings)/mean(df.hh$liq_savings_hat)
df.hh$liq_savings_hat <- df.hh$liq_savings_hat*liqratio

ggplot(df.hh) +
  geom_point(aes(x = log(1+liq_savings), y = liqreg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = log(1+liq_savings), y = log(1+liq_savings_hat), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

#### Housing residuals for depreciation moment
housingfwl <- feols(log(home_value + 1) ~ home_investment| as.factor(hhid) + as.factor(wave),
                 data = df.hh, weights = df.hh$wt_hh, vcov = ~strata)
df.hh$home_value_resid <- housingfwl$residuals
agefwl <- feols(years_ago_built ~ home_investment| as.factor(hhid) + as.factor(wave),
                    data = df.hh, weights = df.hh$wt_hh, vcov = ~strata)
df.hh$years_ago_built_resid <- agefwl$residuals

summary(lm(home_value_resid ~ years_ago_built_resid + 0, data = filter(df.hh, wave!=1), weights = wt_hh))

#### Credit and investment variables adjusted based on ratio of adjusted realized income + expected remittances to actual realized income + expected remittances
df.hh <- mutate(df.hh, E_Y = M_avg + E_remittance_hat,
                credit_hat = credit * E_Y/(total_income + E_remittance_1),
                home_investment_hat = home_investment * E_Y/(total_income + E_remittance_1))

### construct and save adjusted dataset
df.adj <- df.hh %>% 
  select(hhid, wave, wt_hh, food_consumption_hat, home_value_hat, home_investment_hat, liq_savings, tot_savings, E_Y,
         credit, capital_income, credit_cost, credit_hat, liq_savings_hat, tot_savings_hat, total_income_hat, quake_aid, M_avg, max_interest, 
         cut_food, dissaving, E_remittance_hat, avail_remittance_hat, home_value_resid, years_ago_built_resid, default,
         dist_2_seg1pt, quake_losses, received_aid, wt_hh, gorkha_hh, NGO_transfers, designation,
         other_nat_disaster, livestock_farm_shock, riot, illness_injury_shock, price_shock) %>%
  mutate(saved_remittances = (avail_remittance_hat - E_remittance_hat)/E_Y,
         liq_assets = liq_savings_hat/E_Y,
         tot_assets = tot_savings_hat/E_Y,  
         debts = credit_hat/E_Y,
         food_consumption = food_consumption_hat/E_Y,
         home_value = home_value_hat/E_Y,
         home_investment = home_investment_hat/E_Y,
         quake_aid = quake_aid/E_Y,
         total_income = (total_income_hat+E_remittance_hat)/E_Y,
         max_interest = if_else(is.na(max_interest), 0, max_interest)) %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave),
         lag_liq = lag(liq_assets, order_by = wave),
         lag_sav = lag(tot_assets, order_by = wave),
         lag_debt = lag(debts, order_by = wave),
         lag_credit = lag(credit, order_by = wave),
         lag_liq_savings = lag(liq_savings, order_by = wave),
         lag_tot_savings = lag(tot_savings, order_by = wave)) ### raw values for interest rates
saveRDS(df.adj, paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))
write_csv(df.adj, paste(getwd(), "/data/model_output/df_adj.csv", sep = ""))

### Sanity checks
lapply(df.adj, function(x) sum(is.na(x)))
df.adj <- df.adj[complete.cases(df.adj),]

## compute liquidity
df.adj <- mutate(df.adj, liq_test = liq_assets + saved_remittances + food_consumption + home_investment - .05*debts, ## gives beginning of period liquidity including aid
                 liq_test2 = tot_assets + saved_remittances + food_consumption + home_investment - .05*debts)

## compute liquidity
df.adj <- mutate(df.adj, liq_test = 1.05*(lag_liq - lag_debt) + saved_remittances + total_income, ## gives beginning of period liquidity including aid
                 liq_test2 = 1.05*(lag_sav - lag_debt) + saved_remittances + total_income,
                 liq_test3 = 1.05*(lag_liq - lag_debt) + total_income)

ggplot(df.adj, aes(y = log(E_Y), x = ihs(liq_test))) + geom_point(aes(color = log(food_consumption), alpha = .3)) + 
  scale_color_viridis_c() + geom_smooth(method = "lm")
ggplot(df.adj, aes(y = log(E_Y), x = ihs(liq_test2))) + geom_point(aes(color = log(food_consumption), alpha = .3)) + 
  scale_color_viridis_c() + geom_smooth(method = "lm")

### liq savings vs tot savings w/assets
ggplot(df.hh, aes(x = log(1+liq_savings), y = log(1+tot_savings))) + 
  geom_point(alpha = .3, aes(color = gorkha_hh))
cor(df.hh$liq_savings, df.hh$tot_savings, method = "spearman")

ggplot(df.adj) + 
  geom_point(aes(x = ihs(liq_test), y = ihs(lag_h), color = log(food_consumption)), alpha = .3) +
  scale_color_viridis_c()

ggplot(df.adj) + 
  geom_point(aes(x = ihs(liq_test), y = ihs(lag_h), color = home_investment), alpha = .3) +
  scale_color_viridis_c()

ggplot(filter(df.adj, default==1)) + 
  geom_point(aes(x = ihs(liq_test), y = lag_h, color = log(food_consumption))) + scale_color_viridis_c()
ggplot(df.adj) + 
  geom_point(aes(x = ihs(liq_test), y = lag_h, color = implied_default, alpha = .3))

ggplot(df.adj, aes(x = ihs(liq_test), y = food_consumption)) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = ihs(liq_test), y = ihs(home_investment))) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = ihs(lag_h), y = food_consumption)) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = ihs(lag_h), y = ihs(home_investment))) + 
  geom_point() + geom_smooth()

summary(lm(log(food_consumption) ~ ihs(liq_test) + lag_h + ihs(liq_test):lag_h, data = df.adj))
summary(lm(log(1+home_investment) ~ ihs(liq_test) + lag_h + ihs(liq_test):lag_h, data = df.adj))
summary(lm(implied_default ~ ihs(liq_test) + lag_h + ihs(liq_test):lag_h, data = df.adj))

deltmin <- function(theta) {
  df <- df.adj
  delta = theta[1]; R = theta[2]; sigma = theta[3]
  e1 = df$years_ago_built_resid*(df$home_value_resid - (delta-1)*df$years_ago_built_resid)*df$wt_hh/mean(df.hh$wt_hh)
  e2 = ((R-1)*df$liq_savings - df$capital_income)/mean(df$lag_liq_savings)*df$wt_hh/mean(df.hh$wt_hh) 
  e3 = ((R-1)*df$credit - df$credit_cost)/mean(df$lag_credit)*df$wt_hh/mean(df.hh$wt_hh)
  e4 = (sigma^2 * (1 + sigma^2/4) - log(df$total_income)^2)*df$wt_hh/mean(df.hh$wt_hh)
  e = cbind(e1, e2, e3, e4)
  return(sum(colMeans(e)^2))
}

sol <- constrOptim(theta = c(.5, .5, .5), f = deltmin, grad = NULL,
            ui = rbind(diag(3), -1*diag(3)),
            ci = c(rep(0, 3), rep(-5, 3)), 
            control = list(reltol = 3e-2, maxit = 66))
sol
sol <- directL(fn = deltmin,
              lower = c(0,0,0), upper = c(5,5,5), 
              control = list(maxeval = 150))
sol
