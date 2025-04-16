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

df.hh <- read_csv(paste(getwd(), "/data/processed/full_panel.csv", sep = ""), guess_max = 7500) %>% ### drop 1004 observations w missing values
  filter(!is.na(age_hh), !is.na(highest_ed), !is.na(religion), !is.na(home_investment), !is.na(implied_default),
         !is.na(cash_savings), !is.na(loan_payments_recd_ann), !is.na(remittance_income), !is.na(prev_migrants), !is.na(non_durables),
         !is.na(consumption), !is.na(total_income), !is.na(avg_income), !is.na(home_value)) %>%
  mutate(M_avg = NA,
         total_income_hat = NA,
         vdc_id = as.factor(paste(district, vdc, ward, sep = "_")),
         eth_id = as.factor(paste(caste_recode, religion, sep = "_")))

### purge hh heterogeneity and life cycle

avghhdat = data.frame("age_hh" = weighted.mean(df.hh$age_hh, df.hh$wt_hh),
                      "hhmembers" = weighted.mean(df.hh$hhmembers, df.hh$wt_hh),
                      "under18" = weighted.mean(df.hh$under18, df.hh$wt_hh),
                      "over65" = weighted.mean(df.hh$over65, df.hh$wt_hh),
                      "highest_ed" = weighted.mean(df.hh$highest_ed, df.hh$wt_hh),
                      "wave" = 2)

## Expected income regression 
increg <- feols(log(total_income+1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                 poly(highest_ed, 2) + as.factor(wave) +
                 femalehh + poly(landvalue, 2) + prev_migrants + migrant_kat + migrant_india + 
                 migrant_other_asia + migrant_mideast + migrant_earnings | eth_id + vdc_id, 
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
                                                     "migrant_kat" = df.hh$migrant_kat,
                                                     "migrant_india" = df.hh$migrant_india,
                                                     "migrant_other_asia" = df.hh$migrant_other_asia,
                                                     "migrant_mideast" = df.hh$migrant_mideast,
                                                     "migrant_earnings" = df.hh$migrant_earnings,
                                                     "eth_id" = df.hh$eth_id,
                                                     "vdc_id" = df.hh$vdc_id)))
df.hh$total_income_hat <- exp(log(df.hh$M_avg) + increg$residuals)

smearfactor <- weighted.mean(exp(increg$residuals), df.hh$wt_hh)
df.hh$M_avg <- df.hh$M_avg*smearfactor

incratio <- weighted.mean(df.hh$total_income, df.hh$wt_hh)/weighted.mean(df.hh$total_income_hat, df.hh$wt_hh)
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
ggplot() + geom_point(aes(x = log(1+ df.hh$avg_income), y = log(df.hh$M_avg), color = df.hh$gorkha_hh))
## Expected adjusted vals vs realized adjusted vals - how much is permanent vs transitory
ggplot() + geom_point(aes(x = log(df.hh$M_avg), y = log(df.hh$total_income_hat), color = df.hh$gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)


### Purge consumption, housing, default, savings variables
foodreg <- lm(log(non_durables+1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2) + as.factor(wave), 
              data = df.hh, weights = wt_hh)
df.hh$non_durables_hat <- exp(predict(foodreg, newdata = avghhdat) + foodreg$residuals)
foodratio <- weighted.mean(df.hh$non_durables, df.hh$wt_hh)/weighted.mean(df.hh$non_durables_hat, df.hh$wt_hh)
df.hh$non_durables_hat <- df.hh$non_durables_hat*foodratio

summary(foodreg)
ggplot(df.hh) +
  geom_point(aes(x = log(non_durables), y = foodreg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = log(non_durables), y = log(non_durables_hat), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

homereg <- lm(log(home_value+1) ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2) + as.factor(wave), 
              data = df.hh, weights = wt_hh)
df.hh$home_value_hat <- exp(predict(homereg, newdata = avghhdat) + homereg$residuals)
homeratio <- weighted.mean(df.hh$home_value, df.hh$wt_hh)/weighted.mean(df.hh$home_value_hat, df.hh$wt_hh)
df.hh$home_value_hat <- df.hh$home_value_hat*homeratio

summary(homereg)
ggplot(df.hh) +
  geom_point(aes(x = log(home_value + 1), y = homereg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = log(home_value + 1), y = log(home_value_hat), color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)

defreg <- lm(implied_default ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
               poly(highest_ed, 2) + as.factor(wave), 
             data = df.hh, weights = wt_hh)
df.hh$default_hat <- predict(defreg, newdata = avghhdat) + defreg$residuals
defratio <- weighted.mean(df.hh$implied_default, df.hh$wt_hh)/weighted.mean(df.hh$default_hat, df.hh$wt_hh)
df.hh$default_hat <- df.hh$default_hat*defratio

summary(defreg)
ggplot(df.hh) +
  geom_point(aes(x = implied_default, y = defreg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)
ggplot(df.hh) +
  geom_point(aes(x = implied_default, y = default_hat, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0)


#### Housing residuals for depreciation moment
housingfwl <- feols(log(home_value + 1) ~ home_investment| as.factor(hhid) + as.factor(wave),
                    data = df.hh, weights = df.hh$wt_hh, vcov = ~strata)
df.hh$home_value_resid <- housingfwl$residuals
agefwl <- feols(years_ago_built ~ home_investment| as.factor(hhid) + as.factor(wave),
                data = df.hh, weights = df.hh$wt_hh, vcov = ~strata)
df.hh$years_ago_built_resid <- agefwl$residuals

summary(lm(home_value_resid ~ years_ago_built_resid + 0, data = filter(df.hh, wave!=1), weights = wt_hh))

### Liquidity
df.hh$liquidity <- df.hh$lag_cash + df.hh$cap_gains + df.hh$loan_payments_recd_ann + df.hh$remittance_income + df.hh$new_loans_taken + 
  df.hh$inf_transfers + df.hh$total_income + df.hh$quake_aid - df.hh$amount_owed
df.hh$liquidity_w_investments <- df.hh$tot_savings + df.hh$liquidity

hist(df.hh$liquidity[abs(df.hh$liquidity)<2e6], breaks = 100)
hist(df.hh$liquidity_w_investments[abs(df.hh$liquidity_w_investments)<5e6], breaks = 100)
quantile(df.hh$liquidity, na.rm = TRUE)
sum(df.hh$liquidity<0, na.rm = TRUE)
quantile(df.hh$liquidity_w_investments, na.rm = TRUE)
sum(df.hh$liquidity_w_investments<0, na.rm = TRUE)
sum(is.na(df.hh$liquidity) & df.hh$wave!=1) ## new hhs

liqreg <- lm(liquidity ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
               poly(highest_ed, 2) + as.factor(wave), 
             data = df.hh, weights = wt_hh)

df.hh$liquidity_hat <- NA
df.hh$liquidity_hat[!is.na(df.hh$liquidity)] <- predict(liqreg, newdata = avghhdat) + liqreg$residuals
liqratio <- weighted.mean(df.hh$liquidity, df.hh$wt_hh, na.rm = TRUE)/weighted.mean(df.hh$liquidity_hat, df.hh$wt_hh, na.rm = TRUE)
df.hh$liquidity_hat <- df.hh$liquidity_hat*liqratio

summary(liqreg)
ggplot(filter(df.hh, !is.na(liquidity))) + ## 1 outlier
  geom_point(aes(x = liquidity, y = liqreg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0) + xlim(-1e5, 5e6)
ggplot(filter(df.hh, !is.na(liquidity))) +
  geom_point(aes(x = liquidity, y = liquidity_hat, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0) + xlim(-1e5, 1e7) + ylim(-1e5, 1e7)

savereg <- lm(liquidity_w_investments ~ poly(age_hh, 2) + poly(hhmembers, 2) + under18 + over65 +
                poly(highest_ed, 2) + as.factor(wave), 
              data = df.hh, weights = wt_hh)

df.hh$liquidity_plus_hat <- NA
df.hh$liquidity_plus_hat[!is.na(df.hh$liquidity_w_investments)] <- predict(savereg, newdata = avghhdat) + savereg$residuals
savratio <- weighted.mean(df.hh$liquidity_w_investments, df.hh$wt_hh, na.rm = TRUE)/weighted.mean(df.hh$liquidity_plus_hat, df.hh$wt_hh, na.rm = TRUE)
df.hh$liquidity_plus_hat <- df.hh$liquidity_plus_hat*savratio

summary(savereg)
ggplot(filter(df.hh, !is.na(liquidity_w_investments))) +
  geom_point(aes(x = liquidity_w_investments, y = savereg$fitted.values, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0) + xlim(-1e6, 5e7)
ggplot(filter(df.hh, !is.na(liquidity_w_investments))) +
  geom_point(aes(x = liquidity_w_investments, y = liquidity_plus_hat, color = gorkha_hh)) + 
  geom_abline(slope = 1, intercept = 0) + xlim(-1e6, 3e7)

#### Investment variables adjusted based on ratio of adjusted realized income to actual realized income
df.hh$home_investment_hat <- df.hh$home_investment * df.hh$total_income_hat/(df.hh$total_income+1)

df.hh <- mutate(df.hh, total_rainfall = rainfall_month_7 + rainfall_month_8 + rainfall_month_9 + rainfall_month_10 + 
                  rainfall_month_11 + rainfall_month_12 + rainfall_month_1 + rainfall_month_2 + rainfall_month_3 + 
                  rainfall_month_4 + rainfall_month_5 + rainfall_month_6)

### construct and save adjusted dataset
df.adj <- df.hh %>% 
  mutate() %>%
  select(hhid, wave, wt_hh, food_consumption, non_durables_hat, home_value_hat, home_investment_hat, 
         liquidity_plus_hat, liquidity_hat, default_hat, total_income_hat, M_avg, remittance_income,
         years_ago_built_resid, home_value_resid, cash_savings, cap_gains, loan_payments_recd_ann,
         amount_owed, loans_made_1yr, loans_made_2yr, loans_made_3yr, loans_taken_1yr, loans_taken_2yr, loans_taken_3yr, 
         quake_aid, max_interest, cut_food, cut_nonfood, school_interrupt, child_labor, child_wage_labor,
         asset_sales, dissaving, migrants_past_year, dist_2_seg1pt1, dist_2_seg1, quake_losses, gorkha_loss_ever,
         received_aid, wt_hh, gorkha_hh, riot_hh, received_ngo_aid, designation, total_rainfall,
         other_nat_disaster, livestock_farm_shock, riot, illness_injury_shock, price_shock,
         caste_recode, femalehh, landvalue, implied_default, shake_pga, elevation, slope, walls, foundation,
         roof, fuel, stove, currentlyliving, under5, over65) %>%
  mutate(max_interest = if_else(is.na(max_interest), 0, max_interest))

saveRDS(df.adj, paste(getwd(), "/data/model_output/df.adj.rds", sep = ""))

sigmame <- .2
meshocks <- exp(rnorm(nrow(df.adj), sigmame^2/2, sigmame))
df.adj$M_avg <- df.adj$M_avg * meshocks

df.adj <- df.adj %>%
  mutate(liquidity = liquidity_hat/M_avg,
         liquidity_plus = liquidity_plus_hat/M_avg,  
         non_durables = non_durables_hat/M_avg,
         home_value = home_value_hat/M_avg,
         home_investment = home_investment_hat/M_avg,
         quake_aid = quake_aid/M_avg,
         total_income = total_income_hat/M_avg) %>%
  group_by(hhid) %>%
  mutate(lag_h = lag(home_value, order_by = wave),
         lag_y = lag(total_income, order_by = wave))

### Sanity checks
lapply(df.adj, function(x) sum(is.na(x)))
df.adj <- df.adj[complete.cases(df.adj),]

## U shape b/c dividing by low income results in higher liq values
ggplot(df.adj, aes(x = log(M_avg), y = ihs(liquidity))) + geom_point(aes(color = log(non_durables), alpha = .3)) + 
  scale_color_viridis_c() + geom_smooth(method = "lm")
ggplot(df.adj, aes(y = log(M_avg), x = ihs(liquidity_plus))) + geom_point(aes(color = log(non_durables), alpha = .3)) + 
  scale_color_viridis_c() + geom_smooth(method = "lm")

### liq savings vs tot savings w/assets
ggplot(df.hh, aes(x = ihs(liquidity), y = ihs(liquidity_w_investments))) + 
  geom_point(alpha = .3, aes(color = gorkha_hh))
cor(df.hh$liquidity[!is.na(df.hh$liquidity)], df.hh$liquidity_w_investments[!is.na(df.hh$liquidity)], method = "spearman")

ggplot(df.adj) + 
  geom_point(aes(x = ihs(liquidity), y = ihs(lag_h), color = log(non_durables)), alpha = .3) +
  scale_color_viridis_c()

ggplot(df.adj) + 
  geom_point(aes(x = ihs(liquidity), y = ihs(lag_h), color = log(home_investment)), alpha = .3) +
  scale_color_viridis_c()

ggplot(df.adj) + 
  geom_point(aes(x = ihs(liquidity), y = ihs(lag_h), color = default_hat), alpha = .3)

mean(df.adj$non_durables[df.adj$default_hat>.5]); mean(df.adj$non_durables[df.adj$default_hat<.5])
mean(df.adj$home_investment[df.adj$default_hat>.5]); mean(df.adj$home_investment[df.adj$default_hat<.5])
mean(df.adj$liquidity[df.adj$default_hat>.5], na.rm = TRUE); mean(df.adj$liquidity[df.adj$default_hat<.5], na.rm = TRUE)
mean(df.adj$lag_h[df.adj$default_hat>.5], na.rm = TRUE); mean(df.adj$lag_h[df.adj$default_hat<.5], na.rm = TRUE)

ggplot() + geom_density(data = filter(df.adj, default_hat > .5), aes(x = liquidity, color = "red")) +
  geom_density(data = filter(df.adj, default_hat < .5), aes(x = liquidity, color = "blue")) + xlim(-5, 10)

ggplot(df.adj, aes(x = ihs(liquidity), y = non_durables)) +
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = ihs(liquidity), y = ihs(home_investment))) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = ihs(lag_h), y = non_durables)) + 
  geom_point() + geom_smooth()

ggplot(df.adj, aes(x = ihs(lag_h), y = ihs(home_investment))) + 
  geom_point() + geom_smooth()

summary(lm(log(non_durables) ~ ihs(liquidity) + log(lag_h) + ihs(liquidity):log(lag_h), data = df.adj))
summary(lm(log(1+home_investment) ~ ihs(liquidity) + log(lag_h) + ihs(liquidity):log(lag_h), data = df.adj))
summary(lm(default_hat ~ ihs(liquidity) + log(lag_h) + ihs(liquidity):log(lag_h), data = df.adj))

deltmin <- function(theta) {
  df <- df.adj
  delta = theta[1]; R = theta[2]; sigma = theta[3]; sigmame = theta[4]
  e1 = df$years_ago_built_resid*(df$home_value_resid - (delta-1)*df$years_ago_built_resid)*df$wt_hh/mean(df.hh$wt_hh)
  e2 = ((R-1)*df$cash_savings - df$cap_gains)/mean(df$cap_gains)*df$wt_hh/mean(df.hh$wt_hh) 
  e3 = (df$loans_made_1yr*annuitycalc(R, 1) + df$loans_made_2yr*annuitycalc(R, 2) + df$loans_made_3yr*annuitycalc(R, 3) - df$loan_payments_recd_ann)/mean(df$loan_payments_recd_ann)*df$wt_hh/mean(df.hh$wt_hh)
  e4 = (df$loans_taken_1yr*annuitycalc(R, 1) + df$loans_taken_2yr*annuitycalc(R, 2) + df$loans_taken_3yr*annuitycalc(R, 3) - df$amount_owed)/mean(df$amount_owed)*df$wt_hh/mean(df.hh$wt_hh)
  e5 = (log(df$total_income) + sigma^2/2 + sigmame^2/2)*df$wt_hh/mean(df.hh$wt_hh)
  e6 = (log(df$total_income)^2 - (sigma^2/2 + sigmame^2/2)^2 - sigma^2 - sigmame^2)*df$wt_hh/mean(df.hh$wt_hh)
  e = cbind(e1, e2, e3, e4, e5, e6)
  return(sum(colMeans(e)^2))
}

### note this can't really identify sigmame because me shocks are not added within the function
sol <- constrOptim(theta = c(.5, .5, .5, .5), f = deltmin, grad = NULL,
            ui = rbind(diag(4), -1*diag(4)),
            ci = c(rep(0, 4), rep(-5, 4)))
sol
sol <- directL(fn = deltmin,
              lower = c(0,0,0,0), upper = c(5,5,5,5), 
              control = list(maxeval = 2000))
sol

## bound on unpredictable variance of income?
rainreg <- feols(log(total_income) ~ total_rainfall | as.factor(hhid),
                 data = df.adj, weights = df.adj$wt_hh)

sqrt(var(rainreg$fitted.values))
sqrt(var(rainreg$residuals))
