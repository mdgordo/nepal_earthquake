### packages and data
library(tidyverse)
library(parallel)
library(nloptr)
library(fastGHQuad)
library(flatlandr)
library(gmm)
library(Rearrangement)
library(AER)
library(fixest)
options('nloptr.show.inequality.warning'=FALSE)
source(paste(getwd(), "/code/VFIfunctions.R", sep = ""))

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) %>%
  filter(designation !="none", !is.na(age_hh), !is.na(highest_ed), !is.na(inc_wave1), !is.na(home_investment), !is.na(consumption),
         !is.na(credit), !is.na(liq_savings), !is.na(tot_savings), !is.na(remittance_income), !is.na(principal_paid), !is.na(principal_received), 
         !is.na(religion), food_consumption>0, !is.na(total_income), total_income>0) %>%
  mutate(M_avg = NA,
         total_income_hat = NA,
         vdc_id = as.factor(paste(district, vdc, ward, sep = "_")),,
         eth_id = as.factor(paste(caste_recode, religion, sep = "_")))

### purge hh heterogeneity and life cycle
avghhdat = data.frame("age_hh" = mean(df.hh$age_hh),
                      "hhmembers" = mean(df.hh$hhmembers),
                      "under18" = mean(df.hh$under18),
                      "over65" = mean(df.hh$over65),
                      "highest_ed" = mean(df.hh$highest_ed),
                      "wave" = 1)

### income a bit tricky 
increg <- lm(log(total_income) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
               poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh)
df.hh$incresid <- increg$residuals
df.hh$total_income_hat <- exp(predict(increg, avghhdat) + increg$residuals)

## Expected income regression
einc <- feols(incresid ~ femalehh + poly(landvalue, 2) + prev_migrants + migrant_male + migrant_educated +
                migrant_kat + migrant_india + migrant_other_asia + migrant_mideast + migrant_working +
                migrant_earnings | eth_id + vdc_id,
              data = df.hh, weights = df.hh$wt_hh)
df.hh$M_avg <- exp(predict(einc, newdata = data.frame("landvalue" = df.hh$landvalue,
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
                                                     "vdc_id" = df.hh$vdc_id)) + 
                                            predict(increg, avghhdat))

### Expected remittance similarly but tobit
remitreg <- tobit(remittance_income ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                    poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh)
df.hh$remittance_hat <- pmax(0, predict(remitreg, newdata = avghhdat) + resid(remitreg))
df.hh$remittance_resid <- df.hh$remittance_income - pmax(0, predict(remitreg, newdata = df.hh))

erem <- feols(remittance_resid ~ femalehh + poly(landvalue, 2) + prev_migrants + migrant_male + migrant_educated +
                migrant_kat + migrant_india + migrant_other_asia + migrant_mideast + migrant_working +
                migrant_earnings | eth_id + vdc_id,
              data = df.hh, weights = df.hh$wt_hh)
df.hh$E_remittance <- pmax(0, predict(erem, newdata = data.frame("landvalue" = df.hh$landvalue,
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
                                                                 "vdc_id" = df.hh$vdc_id)) + 
                             predict(remitreg, newdata = avghhdat))
df.hh$remittance_loan = pmax(0, df.hh$remittance_hat - df.hh$E_remittance)
df.hh$remittance_saved = pmax(0, df.hh$E_remittance - df.hh$remittance_hat)

### other variables more straightforward

foodreg <- lm(log(food_consumption) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh)
df.hh$food_consumption_hat <- exp(predict(foodreg, newdata = avghhdat) + foodreg$residuals)

## small number of HHs have 0 home value

homereg <- lm(log(home_value+1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                   poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh)
df.hh$home_value_hat <- exp(predict(homereg, newdata = avghhdat) + homereg$residuals)

## tobit for savings, credit, and home investment - larger number of zeros

savereg <- tobit(liq_savings ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                   poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh)
df.hh$liq_savings_hat <- pmax(0, predict(savereg, newdata = avghhdat) + resid(savereg))

savereg2 <- tobit(tot_savings ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                   poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh, control = list(iter.max = 10000))
df.hh$tot_savings_hat <- pmax(0, predict(savereg2, newdata = avghhdat) + resid(savereg2))

credreg <- tobit(credit ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                   poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh)
df.hh$credit_hat <- pmax(0, predict(credreg, newdata = avghhdat) + resid(credreg))

investreg <- tobit(home_investment ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                   poly(highest_ed, 2) + as.factor(wave), data = df.hh, weights = wt_hh)
df.hh$home_investment_hat <- pmax(0, predict(investreg, newdata = avghhdat) + resid(investreg))

df.adj <- df.hh %>% 
  select(hhid, wave, wt_hh, food_consumption_hat, home_value_hat, home_investment_hat, liq_savings, tot_savings, 
         credit, capital_income, credit_cost, credit_hat, liq_savings_hat, tot_savings_hat, total_income_hat, quake_aid, M_avg, max_interest, 
         cut_food, dissaving, E_remittance, remittance_loan, remittance_saved, years_ago_built,
         dist_2_seg1pt, quake_losses, received_aid, wt_hh, gorkha_hh, NGO_transfers, designation,
         other_nat_disaster, livestock_farm_shock, riot, illness_injury_shock, price_shock) %>%
  mutate(E_Y = M_avg + E_remittance,
         liq_assets = liq_savings_hat/E_Y, 
         tot_assets = tot_savings_hat/E_Y,  ## end of period assets
         debts = (credit_hat+remittance_loan)/E_Y,
         food_consumption = food_consumption_hat/E_Y,
         home_value = home_value_hat/E_Y,
         home_investment = home_investment_hat/E_Y,
         quake_aid = quake_aid/E_Y,
         total_income = (total_income_hat+E_remittance)/E_Y,
         years_ago_built = years_ago_built+1) %>%
  group_by(hhid) %>%
  mutate(lag_xp = lag(tot_assets, order_by = wave),
         lag_xn = lag(debts, order_by = wave),
         lag_a = lag(quake_aid, order_by = wave),
         lag_y = lag(total_income, order_by = wave),
         lag_c = lag(food_consumption, order_by = wave),
         lag_h = lag(home_value, order_by = wave),
         lag_i = lag(home_investment, order_by = wave)) 
saveRDS(df.adj, paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))

### Sanity checks

ggplot(df.hh) + geom_point(aes(x = log(M_avg), y = log(total_income)))
ggplot(df.hh) + geom_point(aes(x = log(total_income_hat), y = log(M_avg)))

ggplot(df.hh) + geom_point(aes(x = log(remittance_hat), y = log(remittance_income)))
ggplot(df.hh) + geom_point(aes(x = log(E_remittance), y = log(remittance_income)))
ggplot(df.hh) + geom_histogram(aes(x = remittance_resid), bins = 100)
ggplot(df.hh) + geom_point(aes(x = home_value, y = home_value_hat))
ggplot(df.hh) + geom_point(aes(x = liq_savings, y = liq_savings_hat))
ggplot(df.hh) + geom_point(aes(x = tot_savings, y = tot_savings_hat))
ggplot(df.hh) + geom_point(aes(x = credit, y = credit_hat))
ggplot(df.hh) + geom_point(aes(x = home_investment_hat, y = home_investment))

### liq savings vs tot savings w/assets
ggplot(df.hh, aes(x = log(1+liq_savings), y = log(1+tot_savings))) + 
  geom_point(alpha = .3, aes(color = log(1+gorkha_loss_ever)))
df.hh %>%
  group_by(gorkha_hh) %>%
  summarize(mean(liq_savings),
            mean(tot_savings))
summary(lm(tot_savings ~ liq_savings + log(gorkha_loss_ever+1) + liq_savings:log(gorkha_loss_ever+1), df.hh))

df.adj <- df.adj[complete.cases(df.adj),]
df.adj <- mutate(df.adj, liq_test = tot_assets + food_consumption + home_investment - .05*debts) ## gives beginning of period liquidity including aid
df.adj <- mutate(df.adj, liq_test2 = liq_assets + food_consumption + home_investment - .05*debts)

ggplot(df.adj) + geom_point(aes(x = log(E_Y), y = log(liq_test)))
ggplot(df.adj) + geom_point(aes(x = log(E_Y), y = log(liq_test2)))

ggplot(df.adj) + 
  geom_point(aes(x = liq_test, y = lag_h, color = log(food_consumption)), alpha = .6) +
  scale_color_viridis_c()

ggplot(df.adj) + 
  geom_point(aes(x = liq_test, y = lag_h, color = log(home_investment+1)), alpha = .3) +
  scale_color_viridis_c()

ggplot(df.adj, aes(x = liq_test, y = food_consumption)) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = liq_test, y = home_investment)) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = lag_h, y = food_consumption)) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = lag_h, y = home_investment)) + 
  geom_point() + geom_smooth()

summary(lm(food_consumption ~ liq_test + lag_h + liq_test:lag_h, data = df.adj))
summary(lm(home_investment ~ liq_test + lag_h + liq_test:lag_h, data = df.adj))

deltmin <- function(theta) {
  df <- filter(df.adj, home_value < Inf)
  delta = theta[1]; R = theta[2]; sigma = theta[3]
  e1 = (1/(1-delta) - df$years_ago_built)*df$wt_hh/mean(df.hh$wt_hh)
  e2 = ((R-1)*df$tot_savings - df$capital_income)/mean(df$tot_savings)*df$wt_hh/mean(df.hh$wt_hh) 
  e3 = ((R-1)*df$credit - df$credit_cost)/mean(df$credit)*df$wt_hh/mean(df.hh$wt_hh)
  e4 = (log(df$total_income)^2 - sigma^2)*df$wt_hh/mean(df.hh$wt_hh)
  e = cbind(e1, e2, e3, e4)
  return(sum(colSums(e)^2))
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
