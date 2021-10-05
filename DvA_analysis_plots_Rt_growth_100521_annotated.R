
#############################################################################################################################
# LOAD LIBRARIES & FILES ----------------------------------------------------------------------------------------------------
#############################################################################################################################

library(tidyverse)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(dplyr)
library(EpiEstim)
library(xlsx)
library(zoo)
library(cowplot)
library(data.table)
library(mltools)
library(purrr)
library(DescTools)
library(aweek)
library(devtools)
library(roloc)
library(readr)
library(stats)
library(readxl)
library(readr)
library(covidestim)
library(MASS)

# Load functions

# path <- # path to where functions file is stored
file_name <- "DvA_functions_100521.R"
full_path <- paste0(path, file_name)
source(full_path, echo=TRUE)

# Load GISAID metadata file 
metadata_date <- "2021-08-13"
gisaid_date <- paste("/", metadata_date, "/", sep="")
# path <- # path to where metadata stored 
next_meta <- readNextMeta_states(paste0(path, "metadata_", metadata_date, ".tsv.gz", sep=""))
# save(next_meta, file = "next_meta.Rda") # save so don't have to re-process next time
# load("next_meta.Rda")

# set the number of days for the emergence period
emergence_pd <- 90

#############################################################################################################################
# PREP DATA ----------------------------------------------------------------------------------------------------------------
#############################################################################################################################

# limit to the US
next_meta_sub <- next_meta %>% dplyr::filter(country == "USA")

# limit to states of interest
states <- c("Connecticut", "Massachusetts", "Maine", "New Hampshire", "Rhode Island", "Vermont")
next_meta_sub$state_2 <- trimws(next_meta_sub$state_2) # removes extra spaces before/after state name
next_meta_sub[which(next_meta_sub$state_1 == "CT"), "state_1"] <- "Connecticut"
next_meta_sub[which(next_meta_sub$state_1 == "MA"), "state_1"] <- "Massachusetts"
next_meta_sub[which(next_meta_sub$state_1 == "ME"), "state_1"] <- "Maine"
next_meta_sub[which(next_meta_sub$state_1 == "NH"), "state_1"] <- "New Hampshire"
next_meta_sub[which(next_meta_sub$state_1 == "RI"), "state_1"] <- "Rhode Island"
next_meta_sub[which(next_meta_sub$state_1 == "VT"), "state_1"] <- "Vermont"

# look for state mismatches (state in seqName and Location do not match)
## these states have seqName w/ states of interest but not location
next_meta_sub_seqName <- next_meta_sub %>% dplyr::filter(state_1 %in% states) 

## these states have location w/ states of interest but not seqName
next_meta_sub_location <- next_meta_sub %>% dplyr::filter(state_2 %in% states) 

## not consistent that submitting/originating lab address corresponds to seqName and sample origin corresponds to Location, so dropping those rows
next_meta_sub <- next_meta_sub[which(next_meta_sub$state_1 == next_meta_sub$state_2), ] # also drops any states not in states of interest
next_meta_sub <- next_meta_sub %>% dplyr::filter(state_1 %in% states) # removes XX (probably missing state info)
next_meta_sub$state <- next_meta_sub$state_1
unique(next_meta_sub$state) # check that states look correct
next_meta_NE <- next_meta_sub # keep so can calculate sequencing coverage

# drop sequences w/o lineage assignments
next_meta_sub <- next_meta_sub %>% dplyr::filter(pango_lineage != "None")

# assign to 3 categories
alpha_lineages <- c("B.1.1.7")
delta_lineages <- c("B.1.617.2", paste0("AY.", seq(from = 1, to = 25, by = 1)), "AY.3.1")

next_meta_sub$Variant_Category <- rep(NA, nrow(next_meta_sub))
next_meta_sub[which(next_meta_sub$pango_lineage %in% alpha_lineages), "Variant_Category"] <- "Alpha"
next_meta_sub[which(next_meta_sub$pango_lineage %in% delta_lineages), "Variant_Category"] <- "Delta"
next_meta_sub[which(is.na(next_meta_sub$Variant_Category == TRUE)), "Variant_Category"] <- "Other" 

# check that lineages look correct in each category
unique(next_meta_sub %>% dplyr::filter(Variant_Category == "Alpha") %>% dplyr::select(pango_lineage))
unique(next_meta_sub %>% dplyr::filter(Variant_Category == "Delta") %>% dplyr::select(pango_lineage))
unique(next_meta_sub %>% dplyr::filter(Variant_Category == "Other") %>% dplyr::select(pango_lineage))

# label each variant category w/ binary 0/1
next_states <- encodeLineage_var(next_meta_sub) 

# transform data from wide to long format
next_states_melt <- reshape2::melt(next_states, 
                                   id.vars = c("seqName", "state", "Date"), 
                                   measure.vars =  c("Alpha", "Delta", "Other"))
next_states_melt <- next_states_melt %>% dplyr::rename(Variant_Category = variable)

# select colors for plotting
all_dark2_colors <- brewer.pal(8, "Dark2") 
customPalette <- c(all_dark2_colors[1], all_dark2_colors[3], all_dark2_colors[7]) # Alpha = green; Delta = purple; Other = brown

#############################################################################################################################
# VARIANT FREQUENCIES  -------------------------------------------------------------------------------
#############################################################################################################################

var_categories <- c("Alpha", "Delta", "Other")
freq_time <- plotProp_line(next_states %>% 
                           dplyr::filter(Date >= "2021-01-01", epidate < "2021-08-01"), # restrict dates b/c of delays in sequencing/reporting
                           var_categories)
freq_time

#############################################################################################################################
# LOGISTIC GROWTH DURING VARIANT-SPECIFIC EMERGENCE PERIODS ----------------------------------------------
#############################################################################################################################

# calculate start and end dates of 90 day emergence periods for each variant category / state
var_categories <- c("Alpha", "Delta")
start_end_dates <- calc_date_1st_inf(next_states_melt, var_categories, num_days = emergence_pd)
state_names <- names(start_end_dates)
start_end_dates_rbind <- do.call(rbind.data.frame, start_end_dates)
start_end_dates_rbind$state <- rep(state_names, each = 2)

# restrict data to 90 day emergence periods; add counter
next_states_melt_emerge <- transform_data_Xdays(next_states_melt %>% dplyr::filter(Variant_Category %in% var_categories), 
                                                start_end_dates, 
                                                num_days = emergence_pd)
next_states_melt_emerge <- do.call(rbind.data.frame, next_states_melt_emerge)

# plot binomial logistic growth curves
growth_plots_emerge <- plotLineage_cat_Xdays(next_states_melt_emerge, num_days = emergence_pd)
growth_plots_emerge

# run binomial logistic regression separately so can pull out coefficients
state_biplots <- vector("list", length = length(states))
state_coefs <- vector("list", length = length(states))
predict_list_state <- list()

for (j in 1:length(states)){
  ss_mod <- next_states_melt_emerge %>% dplyr::filter(state == states[j]) 
  mod_size <- ncol(ss_mod)
  ss_mod$Variant_Category <- as.character(ss_mod$Variant_Category)
  model_list <- vector("list", length=length(unique(ss_mod$Variant_Category))) # model output list
  CI_list <- vector("list", length=length(unique(ss_mod$Variant_Category))) # confidence intervals list
  predict_list <- vector("list", length=length(unique(ss_mod$Variant_Category))) # predictions list
  
  for(i in 1:length(model_list)){ # for each variant category
    ss_mod_single <- ss_mod[which(ss_mod$Variant_Category == unique(ss_mod$Variant_Category)[[i]]), ]
    model_list[[i]] <- glm(value ~ Counter, data=ss_mod_single, family = "binomial") # for each state & each lineage, calculate binomial logistic regression
    CI_list[[i]] <- confint(model_list[[i]])[2, 1:2] # calculate confidence intervals around coefficients
    predict_logodds <- predict(model_list[[i]], newdata = data.frame(Counter = c(seq(0, 90, by = 1)), type = "response")) # give fitted values for 1st 90 days
    predict_odds <- exp(predict_logodds)
    predict_probs <- predict_odds / (1 + predict_odds) # transform log odds into probabilities
    predict_list[[i]] <- predict_probs
  }
  
  predict_list_state[[j]] <-  predict_list # fitted values for 1st 90 days
  model_coefs <- sapply(model_list, function(x) coefficients(summary(x))[2,1]) # coefficient for X variable
  coef_std_err <- sapply(model_list, function(x) coefficients(summary(x))[2,2]) # std error for X variable
  model_p <- sapply(model_list, function(x) coefficients(summary(x))[2,4]) # p value for X variable
  CI_vals <- sapply(CI_list, function(x) x[1:2]) # confidence intervals for X variable

    model_df <- data.frame(lineage = unique(ss_mod$Variant_Category),
                         coefs = model_coefs,
                         stderrs = coef_std_err*5.2, # assumes mean generation interval = 5.2
                         pvals = model_p, 
                         logp = -log10(model_p), # takes log10 of p value 
                         transmissibility = exp(model_coefs*5.2),
                         lowci = CI_vals[1, ],
                         upci = CI_vals[2, ])
    
  model_df$logp[is.infinite(model_df$log)] <- 250 
  state_coefs[[j]] <- model_df
}

model_df <- do.call(rbind.data.frame, state_coefs)
model_df$state <- rep(states, each = 2)

# plot multiplicative increase in R 
fold_increase_emerge <- 
  ggplot(model_df, aes(x =transmissibility, y = log(logp), label = lineage)) + 
  geom_point() +
  geom_text(aes(color = lineage), size = 3, position=position_jitter()) + 
  theme_bw() + 
  coord_cartesian(ylim = c(0, 6),  xlim = c(1, 2.2))  +
  labs(x = "Fold increase in R", y = "Evidence (loglogP)") + 
  labs(color = "Variant Category") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  theme(legend.position="none") +
  facet_wrap(~state) +
  scale_y_continuous(breaks = seq(0, 6, by = 1), limits=c(0, 6)) +
  scale_x_continuous(breaks = seq(1, 2.5, by = 0.25), limits=c(1, 2.5))

fold_increase_emerge

# plot binomial logistic regression coefficients bar chart w/ CIs
all_coefs <- do.call(rbind.data.frame, state_coefs)
all_coefs <- all_coefs[, c("lineage", "coefs", "lowci", "upci")]
all_coefs$state <- rep(states, each = length(unique(all_coefs$lineage)))
all_coefs$log_odds <- all_coefs$coefs
all_coefs$lineage <- factor(all_coefs$lineage, levels = c("Delta", "Alpha"))
coef_bar_emerge <- plot_coefs_bar_sep_rotate(all_coefs)
coef_bar_emerge

#############################################################################################################################
#  DAILY VARIANT PROPORTIONS ----------------------------------------------
#############################################################################################################################

# calculate daily proportions of each variant category (note: days w/ no sequences will not be represented)
var_categories <- c("Alpha", "Delta", "Other")
state_vec <- vector("list", length = length(states))
for (j in 1:length(states)){
  pf <- calc_prop_day(next_states %>% dplyr::filter(state == states[j]),
                      var_categories, 
                      "Date") # note: Date = collection date, not GISAID submission date
  state_vec[[j]] <- pf
}
var_data <- state_vec
names(var_data) <- states

# combine into 1 dataframe with states and variant category names
all_state_prop <- list()
for (state in names(var_data)) {
  dat <- var_data[[state]]
  dat_lin_all <- NULL
  for (lineage in names(dat)) {
    dat_lin <- dat[[lineage]]
    dat_lin$lineage <- rep(lineage, nrow(dat_lin))
    dat_lin_all <- rbind.data.frame(dat_lin_all, dat_lin)
  }
  dat_lin_all$state <- rep(state, nrow(dat_lin_all))
  all_state_prop[[state]] <-  dat_lin_all
}
all_state_prop <- do.call(rbind.data.frame, all_state_prop)

# remove Vermont outlier - single sequence after dies out that heavily skews Rt
# View(all_state_prop %>% filter(state == "Vermont", lineage == "Other"))
all_state_prop <- all_state_prop %>% 
  dplyr::mutate(Proportion = ifelse(state == "Vermont" & lineage == "Other" & Date == "2021-07-21", 0, Proportion),
                K = ifelse(state == "Vermont" & lineage == "Other" & Date == "2021-07-21", 0, K),
                N = ifelse(state == "Vermont" & lineage == "Other" & Date == "2021-07-21", 5, N))

all_state_prop <- all_state_prop %>% # have to set N to N-1 for all the other categories; fix proportions in other categories
  dplyr::mutate(Proportion = ifelse(state == "Vermont" & lineage == "Delta" & Date == "2021-07-21", 1, Proportion),
                N = ifelse(state == "Vermont" & lineage %in% c("Alpha", "Delta") & Date == "2021-07-21", 5, N))

# change from long to wide data
all_state_prop_cast_p <- reshape2::dcast(all_state_prop, state + Date ~ lineage, value.var = "Proportion")
all_state_prop_cast_p <- all_state_prop_cast_p %>% dplyr::rename(alpha_prop = Alpha,
                                                                 delta_prop = Delta,
                                                                 other_prop = Other)

# combine daily variant proportions with the number of total sequences (N)
all_state_prop_cast_n <- reshape2::dcast(all_state_prop, state + Date ~ lineage, value.var = "N")
all_state_prop_cast_n <- all_state_prop_cast_n %>%
                         dplyr::select(Date, state, Alpha) %>%  # N same for all lineages, so can just use Alpha
                         dplyr::rename(N = Alpha)
var_data <- left_join(all_state_prop_cast_p, all_state_prop_cast_n, by = c("state" = "state", "Date" = "Date"))

# make sure every date is represented - this is needed for EpiEstim to work correctly
date_range <- rep(seq(as.Date("2021-01-01"), as.Date("2021-08-01"), by = "day"), length(states))
state_range <- rep(states, each = length(date_range)/length(states))
dat_complete_dates <- cbind.data.frame(state_range, date_range)
colnames(dat_complete_dates) <- c("state", "Date")
dat_complete_dates_var <- left_join(dat_complete_dates, var_data, by = c("state" = "state", "Date" = "Date")) 
dat_complete_dates_var[is.na(dat_complete_dates_var)] <- 0 # missing dates (ones w/ no sequences will show as NA - set to 0)
var_data <- dat_complete_dates_var

#############################################################################################################################
# EFFECTIVE REPRODUCTION NUMBER Rt ----------------------------------------------
#############################################################################################################################

# read in and format Covidestim data
# file_name <- # Covidestim csv file
# infect_import = read.csv(file_name)

infect = infect_import %>%
  dplyr::select(state, date, infections, cases.fitted) %>% 
  dplyr::filter(state %in% states, date >= "2021-01-01", date <= "2021-08-01") %>% 
  dplyr::rename(Date = date) %>% 
  dplyr::mutate(Date = as.Date(Date)) 

# calculate 7-day rolling variant proportions
alpha_df_vec <- list()
delta_df_vec <- list()
other_df_vec <- list()

for (state in states) {
  var_data_state <- var_data[which(var_data$state == state), ] 
  var_merge_state <- var_data_state %>% left_join(infect, by = c("state", "Date")) # joins var frequencies w/ estimated infections
  
  var_merge_state <- var_merge_state %>% 
    dplyr::mutate(alpha_infections = infections*alpha_prop, # number of infections per variant category
                  delta_infections = infections*delta_prop,
                  other_infections = infections*other_prop) %>%
    dplyr::mutate(alpha_n = N*alpha_prop, # number of sequences per variant category
                  delta_n = N*delta_prop,
                  other_n = N*other_prop)
  
  var_merge_state <- var_merge_state[order(var_merge_state$Date, decreasing = FALSE), ] # ensures dates in right order for rolling average
  
  daily_7_state <- var_merge_state %>% 
    mutate(alpha_n7 = zoo::rollmean(alpha_n, k = 7, fill = NA), # 7-day rolling avg for # sequences in each variant category
           delta_n7 = zoo::rollmean(delta_n, k = 7, fill = NA),
           other_n7 = zoo::rollmean(other_n, k = 7, fill = NA),
           n_7 = zoo::rollmean(N, k = 7, fill = NA)) %>% # 7-day rolling avg for total # sequences
    mutate(alpha_prop7 = zoo::rollmean(alpha_prop, k = 7, fill = NA), # 7-day rolling avg for variant frequencies (i.e. proportions)
           delta_prop7 = zoo::rollmean(delta_prop, k = 7, fill = NA),
           other_prop7 = zoo::rollmean(other_prop, k = 7, fill = NA)) %>%
    drop_na #drops future dates and first 3 days because of rollmean 
  
  daily_7_state[daily_7_state < 0] <- 0 # set negative values calculated by rolling average to 0 (can't have negative sequences)
  
  alpha_df <- daily_7_state
  delta_df <- daily_7_state
  other_df <- daily_7_state
  
  alpha_df$I <- (alpha_df$alpha_n7 / alpha_df$n_7) * alpha_df$infections # (variant-specific prop) * (# infect) = # variant-specific infections
  delta_df$I <- (delta_df$delta_n7 / delta_df$n_7) * delta_df$infections
  other_df$I <- (other_df$other_n7 / other_df$n_7) * other_df$infections

  # setting NA values to 0; NA comes from n_7 going to 0, then we divide by it above
  alpha_df[is.na(alpha_df)] <- 0
  delta_df[is.na(delta_df)] <- 0
  other_df[is.na(other_df)] <- 0
  
  alpha_df_vec[[state]] <- alpha_df
  delta_df_vec[[state]] <- delta_df
  other_df_vec[[state]] <- other_df
}

# check for outliers - looking for long stretches of zero variant infections 
check_outliers(alpha_df_vec, num_zero_days = 20)
check_outliers(delta_df_vec, num_zero_days = 20)
check_outliers(other_df_vec, num_zero_days = 20)

# Rt method: feeds in time series vector resulting from variant-specific proportion (rolling 7 day avg) * total infections (rolling 7 day avg)
## reflects uncertainty in Rt estimates given the observed seq freq (7 day rolling avg)

alpha_rt_vec <- list()
delta_rt_vec <- list()
other_rt_vec <- list()

MCMC_seed <- 1 # set so can get reproducible results
overall_seed <- 2
mcmc_control <- make_mcmc_control(seed = MCMC_seed)
for (state in states) {
  alpha_rt_vec[[state]] = rt_fun_12inf(alpha_df_vec[[state]], "alpha", window_size = 21, mcmc_control = mcmc_control, seed = overall_seed)
  delta_rt_vec[[state]] = rt_fun_12inf(delta_df_vec[[state]], "delta", window_size = 21, mcmc_control = mcmc_control, seed = overall_seed)
  other_rt_vec[[state]] = rt_fun_12inf(other_df_vec[[state]], "other", window_size = 21, mcmc_control = mcmc_control, seed = overall_seed)
}

# save(alpha_rt_vec, file = "alpha_rt_vec_all.rda") # save if don't want to have to re-run
# save(delta_rt_vec, file = "delta_rt_vec_all.rda")
# save(other_rt_vec, file = "other_rt_vec_all.rda")

# manually truncate Rt estimates plots where start to erroneously increase again despite all zeroes (i.e. variant category died out) due to mean prior of 5

## ALPHA
View(alpha_rt_vec[["Connecticut"]] %>% dplyr::select(Date, rt, I))
alpha_rt_vec[["Connecticut"]] <- alpha_rt_vec[["Connecticut"]] %>% dplyr::filter(Date <= "2021-07-26")

View(alpha_rt_vec[["Maine"]] %>% dplyr::select(Date, rt, I))
alpha_rt_vec[["Maine"]] <- alpha_rt_vec[["Maine"]] %>% dplyr::filter(Date <= "2021-07-26")

View(alpha_rt_vec[["Massachusetts"]] %>% dplyr::select(Date, rt, I))
# alpha_rt_vec[["Massachusetts"]] <- alpha_rt_vec[["Massachusetts"]] %>% dplyr::filter(Date <= "2021-XX-XX")

View(alpha_rt_vec[["New Hampshire"]]%>% dplyr::select(Date, rt, I))
alpha_rt_vec[["New Hampshire"]] <- alpha_rt_vec[["New Hampshire"]] %>% dplyr::filter(Date <= "2021-07-26")

View(alpha_rt_vec[["Rhode Island"]] %>% dplyr::select(Date, rt, I))
alpha_rt_vec[["Rhode Island"]] <- alpha_rt_vec[["Rhode Island"]] %>% dplyr::filter(Date <= "2021-07-23")

View(alpha_rt_vec[["Vermont"]] %>% dplyr::select(Date, rt, I))
alpha_rt_vec[["Vermont"]] <- alpha_rt_vec[["Vermont"]] %>% dplyr::filter(Date <= "2021-07-22")

## DELTA
# View(delta_rt_vec[["Connecticut"]] %>% dplyr::select(Date, rt, I))
# View(delta_rt_vec[["Maine"]] %>% dplyr::select(Date, rt, I))
# View(delta_rt_vec[["Massachusetts"]] %>% dplyr::select(Date, rt, I))
# View(delta_rt_vec[["New Hampshire"]] %>% dplyr::select(Date, rt, I))
# View(delta_rt_vec[["Rhode Island"]] %>% dplyr::select(Date, rt, I))
# View(delta_rt_vec[["Vermont"]] %>% dplyr::select(Date, rt, I))

## OTHER
View(other_rt_vec[["Connecticut"]] %>% dplyr::select(Date, rt, I))
# other_rt_vec[["Connecticut"]] <- other_rt_vec[["Connecticut"]] %>% dplyr::filter(Date <= "2021-XX-XX")

View(other_rt_vec[["Maine"]] %>% dplyr::select(Date, rt, I))
# other_rt_vec[["Maine"]] <- other_rt_vec[["Maine"]] %>% dplyr::filter(Date <= "2021-XX-XX")

View(other_rt_vec[["Massachusetts"]] %>% dplyr::select(Date, rt, I))
# other_rt_vec[["Massachusetts"]] <- other_rt_vec[["Massachusetts"]] %>% dplyr::filter(Date <= "2021-XX-XX")

View(other_rt_vec[["New Hampshire"]] %>% dplyr::select(Date, rt, I))
other_rt_vec[["New Hampshire"]] <- other_rt_vec[["New Hampshire"]] %>% dplyr::filter(Date <= "2021-07-27")

View(other_rt_vec[["Rhode Island"]] %>% dplyr::select(Date, rt, I))
# other_rt_vec[["Rhode Island"]] <- other_rt_vec[["Rhode Island"]] %>% dplyr::filter(Date <= "2021-XX-XX")

View(other_rt_vec[["Vermont"]] %>% dplyr::select(Date, rt, I))
other_rt_vec[["Vermont"]] <- other_rt_vec[["Vermont"]] %>% dplyr::filter(Date <= "2021-07-01")

# binda nd save truncated results together
alpha_rt <- add_state_rbind(alpha_rt_vec, states)
delta_rt <- add_state_rbind(delta_rt_vec, states)
other_rt <- add_state_rbind(other_rt_vec, states)

# save(alpha_rt, file = "alpha_rt_I_trunc.rda") # save truncated data for plotting
# save(delta_rt, file = "delta_rt_I_trunc.rda")
# save(other_rt, file = "other_rt_I_trunc.rda")

# plot Rt estimates against infections to visually inspect for irregular behavior
## ALPHA
scale_y <- mean(alpha_rt$I / alpha_rt$rt, na.rm = TRUE)
ggplot(alpha_rt) + geom_line(aes(Date, I), color = "red") + geom_line(aes(x=Date, y = rt * scale_y)) +
  scale_y_continuous(name = "Infections", sec.axis = sec_axis(~./scale_y, name = "rt")) + facet_wrap(~state) +
  theme(axis.title.y = element_text(color = "red"),
        axis.title.y.right = element_text(color = "black")) 

ggplot(alpha_rt) + geom_line(aes(x=Date, y = rt)) + facet_wrap(~state) +
  geom_hline(yintercept = 1, linetype = 1, color = "red", size = 0.3)

ggplot(alpha_rt) + geom_line(aes(x=Date, y = I)) + facet_wrap(~state, scales = "free_y")

## DELTA
scale_y <- mean(delta_rt$I / delta_rt$rt, na.rm = TRUE)
ggplot(delta_rt) + geom_line(aes(Date, I), color = "red") + geom_line(aes(x=Date, y = rt * scale_y)) +
  scale_y_continuous(name = "Infections", sec.axis = sec_axis(~./scale_y, name = "rt")) + facet_wrap(~state) +
  theme(axis.title.y = element_text(color = "red"),
        axis.title.y.right = element_text(color = "black"))

ggplot(delta_rt) + geom_line(aes(x=Date, y = rt)) + facet_wrap(~state) +
  geom_hline(yintercept = 1, linetype = 1, color = "red", size = 0.3)

ggplot(delta_rt) + geom_line(aes(x=Date, y = I)) + facet_wrap(~state, scales = "free_y")

## OTHER
scale_y <- mean(other_rt$I / other_rt$rt, na.rm = TRUE)
ggplot(other_rt) + geom_line(aes(Date, I), color = "red") + geom_line(aes(x=Date, y = rt * scale_y)) +
  scale_y_continuous(name = "Infections", sec.axis = sec_axis(~./scale_y, name = "rt")) + facet_wrap(~state) +
  theme(axis.title.y = element_text(color = "red"),
        axis.title.y.right = element_text(color = "black"))

ggplot(other_rt) + geom_line(aes(x=Date, y = rt)) + facet_wrap(~state) +
  geom_hline(yintercept = 1, linetype = 1, color = "red", size = 0.3)

ggplot(other_rt) + geom_line(aes(x=Date, y = I)) + facet_wrap(~state, scales = "free_y")

# prepare to plot Rt 
alpha_rt$Variant_Category <- rep("Alpha", nrow(alpha_rt))
delta_rt$Variant_Category <- rep("Delta", nrow(delta_rt))
other_rt$Variant_Category <- rep("Other", nrow(other_rt))

alpha_rt <- alpha_rt %>% 
  dplyr::select(Date, state, Variant_Category, rt, rtlowci, rtupci)

delta_rt <- delta_rt %>% 
  dplyr::select(Date, state, Variant_Category, rt, rtlowci, rtupci)

other_rt <- other_rt %>% 
  dplyr::select(Date, state, Variant_Category, rt, rtlowci, rtupci)

all_rt <- rbind.data.frame(alpha_rt, delta_rt, other_rt)
                      
all_rt$Variant_Category <- factor(all_rt$Variant_Category, levels = c("Alpha", "Delta", "Other"))

# check where to set limits for plotting
min_plot <- min(all_rt$rtlowci, na.rm = TRUE)
max_plot_upci <- max(all_rt$rtupci, na.rm = TRUE)
max_plot_mean <- max(all_rt$rt, na.rm = TRUE)

# chose to restrict upper limit so report those in paper
maine_rt <- all_rt %>% filter(state == "Maine")
round(max(maine_rt$rtupci, na.rm = TRUE), digits = 2)

vermont_rt <- all_rt %>% filter(state == "Vermont")
round(max(vermont_rt$rtupci, na.rm = TRUE), digits = 2)

rt_plot <- ggplot(all_rt, aes(Date, rt)) +
  geom_line(aes(color = Variant_Category)) +
  geom_ribbon(aes(ymin=rtlowci, ymax=rtupci, fill = Variant_Category), linetype = 0, alpha = 0.2) +
  scale_fill_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) +
  scale_color_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) + 
  scale_y_continuous(breaks = seq(0, max_plot_mean, by = 1), limits = c(0, max_plot_mean)) +
  facet_wrap(~state, ncol = 1) +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Rt") +
  theme(legend.position = "bottom", legend.box = "horizontal")  +
  labs(color = "Variant Category", fill = "Variant Category") 

rt_plot

#############################################################################################################################
# EFFECTIVE REPRODUCTION NUMBER Rt RATIOS ----------------------------------------------
#############################################################################################################################

# calculate Rt ratio - Delta:Alpha 
delta_to_alpha <- left_join((all_rt %>% dplyr::filter(Variant_Category == "Delta")), 
                            (all_rt %>% dplyr::filter(Variant_Category == "Alpha")),
                            by = c("state", "Date"))
delta_to_alpha$Rt_Ratio <- delta_to_alpha$rt.x / delta_to_alpha$rt.y # Delta Rt / Alpha Rt
delta_to_alpha$Comparison <- "Delta:Alpha"

# calculate Rt ratio - Delta:Other
delta_to_other <- left_join((all_rt %>% dplyr::filter(Variant_Category == "Delta")), 
                            (all_rt %>% dplyr::filter(Variant_Category == "Other")),
                            by = c("state", "Date"))
delta_to_other$Rt_Ratio <- delta_to_other$rt.x / delta_to_other$rt.y
delta_to_other$Comparison <- "Delta:Other"

delta_comp_all <- rbind.data.frame(delta_to_alpha, delta_to_other)

# start at the earliest date there is a Delta variant + 21 days (b/c of 21 window for estimation) - o/w bunch of NA
first_delta_date <- min(start_end_dates_rbind[which(start_end_dates_rbind$var == "Delta"), "min_date_var"]) + 21
delta_comp_all_limit <- delta_comp_all %>% dplyr::filter(Date >= first_delta_date)

# set as factor for plotting
delta_comp_all_limit$Comparison <- factor(delta_comp_all_limit$Comparison, levels = c("Delta:Alpha", "Delta:Other"))

# plot Rt ratios - Delta:Alpha only
delta_comp_all_limit_ad <- delta_comp_all_limit %>% dplyr::filter(Comparison == "Delta:Alpha")
max_level <- max(delta_comp_all_limit_ad$Rt_Ratio, na.rm = TRUE)

rt_ratio_ad <- delta_comp_all_limit_ad %>%
  ggplot(aes(x = Date, y = Rt_Ratio, color = Comparison, fill = Comparison)) +
  geom_line() +
  scale_fill_manual(values = "dark grey") +
  scale_color_manual(values = "dark grey") +
  facet_wrap(~state, ncol = 1) +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Rt Ratio")+
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_y_continuous(breaks = seq(0, max_level, by = 1), limits = c(0, max_level)) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme(legend.position = "bottom", legend.box = "horizontal") 

rt_ratio_ad

#############################################################################################################################
# DOTPLOT Ratios ----------------------------------------------
#############################################################################################################################

#  plotting all daily Rt ratios for each comparison - Delta:Alpha only
delta_comp_all_limit_ad <- delta_comp_all_limit %>% dplyr::filter(Comparison == "Delta:Alpha")

rt_ratio_dotplot_ad <- 
  ggplot(delta_comp_all_limit_ad, aes(x=Comparison, y=Rt_Ratio, color = Comparison)) +
  geom_point(aes(color = Comparison), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = "dark grey") +  
  scale_fill_manual(values = "dark grey", guide = guide_legend(override.aes = list(alpha = 0))) +
  xlab("Comparison - Delta:Alpha") +
  ylab("Daily Estimated Rt Ratios") +
  theme_bw() +
  facet_wrap(~state, ncol = 1) +
  scale_y_continuous(breaks = seq(0, 4.7, by = 1), limits = c(0, 4.7)) +
  theme(legend.position = "bottom", 
        legend.box = "horizontal",
        legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent"),
        axis.text.y = element_blank()) +
  guides(colour=guide_legend(override.aes=list(size=0, color = "transparent"))) +
  coord_flip()

rt_ratio_dotplot_ad

#############################################################################################################################
# Percent of Confirmed Cases Sequenced over Time ----------------------------------------------
#############################################################################################################################

library(dplyr)
library(plyr)

# read in JHU confirmed cases
urlfile="https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv"
dat_cases <- read_csv(url(urlfile))

dat_cases <- dat_cases %>%
  dplyr::filter(Province_State %in% states) %>%
  dplyr::select(-c(iso2, iso3, code3, FIPS, Admin2, Country_Region, Lat, Long_, Combined_Key)) # drop columns

dat_cases_melt <- reshape2::melt(dat_cases, id.vars = c("UID", "Province_State"))
dat_cases_melt_sum <- ddply(dat_cases_melt, .(Province_State, variable), summarize, Cases_Sum = sum(value, na.rm = TRUE)) # variable = date
dat_cases_melt_sum <- dat_cases_melt_sum %>%
  dplyr::rename(state = Province_State,
                Date = variable,
                Cases = Cases_Sum)
dat_cases_melt_sum$Date <- as.Date(dat_cases_melt_sum$Date, format = "%m/%d/%y")
cases <- dat_cases_melt_sum

# calculate daily incidence from cumulative cases
library(sars2pack)
cases_daily <- add_incidence_column(cases, date_column = "Date", count_column = "Cases", incidence_col_name = "Cases_Daily", grouping_columns = "state")
detach("package:sars2pack", unload = TRUE)

# calculate number of daily sequences
seq <- ddply(next_meta_NE, .(state, Date), nrow) # number of seq per state by date; note includes "None" assignments - just to get a sense of how much seq being done
seq <- seq %>% dplyr::rename(Num_Genomes = V1)

# combine incident daily cases and daily sequences
cases_seq <- left_join(cases_daily, seq, by = c("state" = "state", "Date" = "Date"))
cases_seq_date <- cases_seq %>%
  dplyr::filter(Date >= "2021-01-01") %>%
  dplyr::rename(Genomes_Daily = Num_Genomes)
cases_seq_date[is.na(cases_seq_date)] <- 0 # for days where no sequences

cases_seq_list <- list()
for (state in states) {
  cases_seq_date_state <- cases_seq_date[which(cases_seq_date$state == state), ]
  cases_seq_date_state <- cases_seq_date_state[order(cases_seq_date_state$Date, decreasing = FALSE), ]
  cases_seq_date_state$Cases_Daily7 <- zoo::rollmean(cases_seq_date_state$Cases_Daily, k = 7, fill = NA)
  cases_seq_date_state$Genomes_Daily7 <- zoo::rollmean(cases_seq_date_state$Genomes_Daily, k = 7, fill = NA)
  cases_seq_date_state$Percent_Seq7 = cases_seq_date_state$Genomes_Daily7/cases_seq_date_state$Cases_Daily7
  cases_seq_list[[state]] <- cases_seq_date_state
}

cases_seq_date_all <- do.call(rbind.data.frame, cases_seq_list)

# add 2019 state populations
dat_pop <- read_excel("nst-est2019-01.xlsx")
dat_pop <- dat_pop[4:nrow(dat_pop), c(1, 13)]
colnames(dat_pop) <- c("state", "pop")
dat_pop$state <- gsub("[.]", "", dat_pop$state)
dat_pop <- dat_pop %>% dplyr::filter(state %in% states)

cases_seq_date_7 <- cases_seq_date_all %>% dplyr::left_join(dat_pop, by = "state") 

cases_seq_list7 <- list()
for (state in states) {
  cases_seq_date_7_state <- cases_seq_date_7[which(cases_seq_date_7$state == state), ]
  cases_seq_date_7_state$Cases_Daily_per_100K7 <- (cases_seq_date_7_state$Cases_Daily7 / cases_seq_date_7_state$pop) * 100000
  cases_seq_list7[[state]] <- cases_seq_date_7_state
}
cases_seq_date_7_all <- do.call(rbind.data.frame, cases_seq_list7)

# plot percent of cases sequenced
min_date <- as.Date("2021-01-01")
max_date <- as.Date("2021-08-01")

sec_axis_scale <- max(cases_seq_date_7_all$Cases_Daily_per_100K7, na.rm = TRUE) / max(cases_seq_date_7_all$Percent_Seq7, na.rm = TRUE)

perc_cases_seq_roll7 <- 
  ggplot(cases_seq_date_7_all) + 
  geom_bar(aes(x = Date, y = Cases_Daily_per_100K7, fill = state), stat = "identity") +
  geom_line(aes(x = Date, y = Percent_Seq7*sec_axis_scale, color = state)) +
  scale_color_manual(values = rep("black", 6)) +
  scale_fill_manual(values = rep("grey", 6)) +
  theme_bw() +
  scale_y_continuous(name = "Cases per 100K Population", 
                     sec.axis = sec_axis(~./sec_axis_scale, name = "Percent of Cases Sequenced", 
                     labels = function(b) { paste0(round(b * 100, 0), "%")})) + 
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")) + 
  facet_wrap(~state, ncol = 1) +
  xlim(min_date, max_date) +  
  theme(legend.position = "bottom", legend.box = "horizontal") +
  theme(legend.title = element_text(color = "transparent"),
        legend.text = element_text(color = "transparent")) 
  
perc_cases_seq_roll7

#############################################################################################################################
# Number of Genomes by State, Variant Category ----------------------------------------------
#############################################################################################################################

dat <- data.frame(cases_seq_date_7_all) %>% dplyr::filter(Date <= as.Date("2021-08-01"))
ggplot(dat, aes(Date, Genomes_Daily7)) + geom_line() + facet_wrap(~state)
ggplot(dat, aes(Date, Cases_Daily7)) + geom_line() + facet_wrap(~state)

next_states_melt_sum <- ddply(next_states_melt, .(state, Variant_Category, Date), summarize, Num_Genomes = sum(value, na.rm = TRUE))
next_states_melt_sum$Num_Genomes_7 <- zoo::rollmean(next_states_melt_sum$Num_Genomes, k = 7, fill = NA)
dat <- next_states_melt_sum %>% dplyr::filter(Date >= "2021-01-01", Date <= "2021-08-01")
ggplot(dat, aes(Date, Num_Genomes_7, color = Variant_Category)) + geom_line() + facet_wrap(~state, scales = "free_y")

#############################################################################################################################
# Percent Fully Vaccinated Over Time ----------------------------------------------
#############################################################################################################################

# get vx prop by state over time
# https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc 
dat_vx <- read_csv("COVID-19_Vaccinations_in_the_United_States_Jurisdiction.csv")

dat_vx[which(dat_vx$Location == "MA"), "Location"] <- "Massachusetts"
dat_vx[which(dat_vx$Location == "RI"), "Location"] <- "Rhode Island"
dat_vx[which(dat_vx$Location == "CT"), "Location"] <- "Connecticut"
dat_vx[which(dat_vx$Location == "ME"), "Location"] <- "Maine"
dat_vx[which(dat_vx$Location == "NH"), "Location"] <- "New Hampshire"
dat_vx[which(dat_vx$Location == "VT"), "Location"] <- "Vermont"

dat_vx_sub <- dat_vx %>% 
  dplyr::select(Date, MMWR_week, Location, Series_Complete_Pop_Pct) %>%
  dplyr::filter(Location %in% states) %>%
  dplyr::rename(state = Location) %>%
  dplyr::mutate(Date = as.Date(Date, format = '%m/%d/%Y'))

dat_vx_sub <- dat_vx_sub %>% dplyr::filter(Date >= "2021-01-01")
ggplot(dat_vx_sub, aes(Date, Series_Complete_Pop_Pct, color = state)) + geom_line()

# combine vx w/ confirmed infections over time
dat_vx_sub_cases <- dat_vx_sub %>% dplyr::left_join(infect, by = c( "state", "Date"))

dat_vx_state_list <- list()
for (state in states) {
  dat_vx_case_state <- dat_vx_sub_cases[which(dat_vx_sub_cases$state == state), ]
  dat_vx_case_state <- dat_vx_case_state[order(dat_vx_case_state$Date, decreasing = FALSE), ]
  dat_vx_case_state$Series_Complete_Pop_Pct_7 <- zoo::rollmean(dat_vx_case_state$Series_Complete_Pop_Pct, k = 7, fill = NA)
  dat_vx_state_list[[state]] <- dat_vx_case_state
}

dat_vx_cases_all <- do.call(rbind.data.frame, dat_vx_state_list)
dat_vx_cases_all <- dat_vx_cases_all %>% left_join(dat_pop, by = c("state"))
dat_vx_cases_all$Inf_per_100K <- (dat_vx_cases_all$infections / dat_vx_cases_all$pop) * 100000

# plot perc vx over time & # of cases
vx_inf <- plot_vx_inf(dat_vx_cases_all, start_end_dates_rbind) 
vx_inf_labels <- data.frame(state = c("Connecticut"), 
                            label = c("Alpha Emergence Period"))
vx_inf <- vx_inf + geom_text(x = as.Date("2021-02-15"), y = 275, color = "black", size = 2.5, aes(label = label), data = vx_inf_labels)
vx_inf_labels <- data.frame(state = c("Connecticut"), 
                            label = c("Delta Emergence Period"))
vx_inf <- vx_inf + geom_text(x = as.Date("2021-06-03"), y = 275, color = "black", size = 2.5, aes(label = label), data = vx_inf_labels)
vx_inf

#############################################################################################################################
# Figure 1: Sequencing Coverage & Variant Frequencies ----------------------------------------------
#############################################################################################################################

p <- ggarrange(perc_cases_seq_roll7 + guides(color = guide_legend(nrow=2,byrow=TRUE,override.aes = list(alpha = 0)),
                                             fill=guide_legend(nrow=2,byrow=TRUE,override.aes = list(alpha = 0))),
               freq_time + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
               labels = c("A", "B"),
               nrow = 1,
               widths = c(1.5, 1.5), 
               heights = c(5,5))

ggsave(paste(path, "figures", gisaid_date, "Fig_1.jpg", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Figure 2: Vx/Inf, Growth Rates, Slope, Vx & Inf vs Rel Growth Rates ----------------------------------------------
#############################################################################################################################

library(gridExtra)
library(ggpubr)

mylegend <- g_legend(growth_plots_emerge)
p <- grid.arrange(ggarrange(vx_inf + theme(legend.position="none"),
                            growth_plots_emerge + theme(legend.position="none"),
                            coef_bar_emerge + theme(legend.position="none"),
                            nrow = 1,
                            widths = c(1.3, 0.85, 0.85),
                            labels = c("A", "B", "C")),
                  mylegend,
                  nrow = 2,
                  heights = c(9.5, 0.5))

ggsave(paste(path, "figures", gisaid_date, "Fig_2.jpg", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Figure 3: Rt Line Plots, Ratios, Dotplots  ----------------------------------------------
#############################################################################################################################

p <- ggarrange(rt_plot +guides(fill=guide_legend(nrow=2,byrow=TRUE)), 
               rt_ratio_ad + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
               rt_ratio_dotplot_ad + guides(color=guide_legend(nrow=2,byrow=TRUE, override.aes = list(alpha = 0, size = 0))),
               labels = c("A", "B", "C"),
               ncol = 3,
               nrow = 1,
               widths = c(1.5, 1.5, 1.1), 
               heights = c(5, 5, 5))

ggsave(paste(path, "figures", gisaid_date, "Fig_3.jpg", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Supplementary Figure 1A-C: Vx/Inf/Date of Detection  ----------------------------------------------
#############################################################################################################################

# get vx and inf numbers for the date of 1st delta emergence in each state
delta_start <- start_end_dates_rbind %>% 
  dplyr::filter(var == "Delta") %>% 
  dplyr::select(state, min_date_var) %>% 
  dplyr::rename(Date = min_date_var)
dat_delta_start_vx_inf <- delta_start %>% dplyr::left_join(dat_vx_cases_all, by = c("state", "Date"))

# add relative logistic regression coefficients: Delta / Alpha
all_coefs_cast <- reshape2::dcast(all_coefs, state ~ lineage, value.var = "log_odds")
all_coefs_cast$delta_div_alpha_perc <- ((all_coefs_cast$Delta - all_coefs_cast$Alpha)/all_coefs_cast$Alpha) * 100
all_coefs_cast$delta_div_alpha <- (all_coefs_cast$Delta/all_coefs_cast$Alpha)

dat_delta_start_vx_inf_coef <- dat_delta_start_vx_inf %>% dplyr::left_join(all_coefs_cast, by = "state")

emerge_gr_vx <- ggplot(dat_delta_start_vx_inf_coef, aes(x = Series_Complete_Pop_Pct_7, y = delta_div_alpha, label = state)) +
  geom_point() +
  geom_text(size = 3, nudge_x = 1.4, nudge_y = -0.05) +
  xlab("Percent of Population Fully Vaccinated at 1st Delta Detection") +
  ylab("Logistic Growth Rate - Delta:Alpha")

emerge_gr_inf <- ggplot(dat_delta_start_vx_inf_coef, aes(x = Inf_per_100K, y = delta_div_alpha, label = state)) +
  geom_point() +
  geom_text(size = 3, nudge_x = 2, nudge_y = -0.05) +
  xlab("Infections per 100K Population at Time of 1st Delta Detection") +
  ylab("Logistic Growth Rate - Delta:Alpha")

emerge_vx_date <- ggplot(dat_delta_start_vx_inf_coef, aes(x = Date, y = Series_Complete_Pop_Pct, label = state)) +
  geom_point() +
  geom_text(size = 3, nudge_x = 1, nudge_y = -0.5) +
  xlab("Date of 1st Delta Detection") +
  ylab("Percent of Population Fully Vaccinated at 1st Delta Detection")

p <- grid.arrange(ggarrange(emerge_gr_vx + theme(legend.position="none"),
                            emerge_gr_inf + theme(legend.position="none") ,
                            emerge_vx_date + theme(legend.position="none"),
                            labels = c("A", "B", "C")))

ggsave(paste(path, "figures", gisaid_date, "Fig_S1.jpg", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Supplementary Figure 2: Multiplicative Increase in R ----------------------------------------------
#############################################################################################################################

p <- fold_increase_emerge
ggsave(paste(path, "figures", gisaid_date, "Fig_S2.jpg", sep=""), p, height = 5, width =10)

#############################################################################################################################
# Summary Stats  ----------------------------------------------
#############################################################################################################################

# Sequencing Coverage ***********************************************************************

## maximum weekly coverage by state
max_seq_cov_state <- ddply(cases_seq_date_7_all, .(state), summarize, Max_Seq = round(max(Percent_Seq7, na.rm = TRUE), digits = 2))
max_seq_cov_state[order(max_seq_cov_state$Max_Seq, decreasing = TRUE), ]
mean_cov_state <- ddply(cases_seq_date_7_all, .(state), summarize, Mean_Seq = mean(Percent_Seq7, na.rm = TRUE))

## date all states hit at least 5% seq coverage consistently
above_5p_seq_cov <- cases_seq_date_7_all
above_5p_seq_cov$Above_5 <- rep(NA, nrow(above_5p_seq_cov))
above_5p_seq_cov[which(above_5p_seq_cov$Percent_Seq7 >= 0.05), "Above_5"] <- "Y"
above_5p_seq_cov[which(above_5p_seq_cov$Percent_Seq7 < 0.05), "Above_5"] <- "N"

above_5p_seq_cov_first_time <- ddply(above_5p_seq_cov %>% dplyr::filter(Above_5 == "Y"), .(state), summarize, Min_5perc = min(Date, na.rm = TRUE))
above_5p_seq_cov <- above_5p_seq_cov %>% filter(Date >= max(above_5p_seq_cov_first_time$Min_5perc))
any(above_5p_seq_cov$Percent_Seq7 < 0.05) # check that no states dipped back below 0.05
check_5perc <- above_5p_seq_cov %>% dplyr::filter(Above_5 == "N", Date < "2021-07-15")
ddply(check_5perc, .(state), summarize, Max_Date_Below_5p = max(Date))

# Variant Frequencies ***********************************************************************
df <- NULL
for (lineage in var_categories) {
  l_int <- next_states %>%
    dplyr::group_by(state, epidate) %>%
    dplyr::summarise(Proportion = mean(.data[[lineage]], na.rm=T), K= sum(.data[[lineage]], na.rm=T), N = n())
  l_int$variant_category <- rep(lineage, nrow(l_int))
  df <- rbind.data.frame(df, l_int)
}

df <- as.data.frame(df)

## alpha freq
df_alpha <- df %>% dplyr::filter(variant_category == "Alpha")
max_freq_alpha <- ddply(df_alpha, .(state), summarize, Max_Alpha_Freq = round(max(Proportion), digits = 2))
max_freq_alpha[order(max_freq_alpha$Max_Alpha_Freq, decreasing = TRUE), ]

## delta freq
df_delta <- df %>% dplyr::filter(variant_category == "Delta", epidate == "2021-07-26")
df_delta$Proportion <- round(df_delta$Proportion, digits = 2)
df_delta[order(df_delta$Proportion, decreasing = TRUE), ]

# Logistic Growth Rates  ***********************************************************************

## prob of a given sequence belonging to Alpha v Delta at the start of emergence period
## date reached >50% predicted freq
## number days to dominance (Delta)
dat_full <- NULL
for (state in 1:length(predict_list_state)) {
  
  dat_state <- predict_list_state[[state]]
  dat_var <- NULL
  
  for (variant in 1:length(dat_state)) {
    dat <- NULL
    dat_state_var <- dat_state[[variant]]
    dat$State <- states[state]
    
    if (variant == 1) {
      dat$Variant_Category <- "Alpha" } else {
        dat$Variant_Category <- "Delta"
      }
    dat$Starting_Value <- dat_state_var[1]
    dom_value <- max(which(dat_state_var < 0.5)) + 1
    dat$Over_50_Value <- dom_value
    dat_var <- rbind.data.frame(dat_var, dat)
  }
  dat_full <- rbind.data.frame(dat_full, dat_var)
  dat_full
}

dat_full$Starting_Value <- round(dat_full$Starting_Value, digits = 2)
dat_full_starting_val <- reshape2::dcast(dat_full, State ~ Variant_Category, value.var = "Starting_Value")
dat_full_starting_val[order(dat_full_starting_val$Alpha, decreasing = TRUE), ]

dat_full_dom_val <- reshape2::dcast(dat_full, State ~ Variant_Category, value.var = "Over_50_Value")
dat_full_dom_val[order(dat_full_dom_val$Delta, decreasing = TRUE), ]
round(mean(dat_full_dom_val$Delta), digits = 0)

## relative growth rates: Delta to Alpha
all_coefs_cast$delta_div_alpha_perc <- round(all_coefs_cast$delta_div_alpha_perc, digits = 0)
all_coefs_cast$delta_div_alpha <- round(all_coefs_cast$delta_div_alpha, digits = 2)

all_coefs_cast[order(all_coefs_cast$delta_div_alpha_perc, decreasing = TRUE), ]

# Vaccination % Complete at Start of Alpha vs Delta  ***********************************************************************
start_dates_alpha <- start_end_dates_rbind %>% dplyr::filter(var == "Alpha") 
start_dates_alpha_vx <- left_join(start_dates_alpha, dat_vx_cases_all, by = c("state" = "state", "min_date_var" = "Date"))
start_dates_alpha_vx$Series_Complete_Pop_Pct_7 <- round(start_dates_alpha_vx$Series_Complete_Pop_Pct_7, digits = 2)
start_dates_alpha_vx[order(start_dates_alpha_vx$Series_Complete_Pop_Pct_7, decreasing = TRUE), ]

start_dates_delta <- start_end_dates_rbind %>% dplyr::filter(var == "Delta") 
start_dates_delta_vx <- left_join(start_dates_delta, dat_vx_cases_all, by = c("state" = "state", "min_date_var" = "Date"))
start_dates_delta_vx$Series_Complete_Pop_Pct_7 <- round(start_dates_delta_vx$Series_Complete_Pop_Pct_7, digits = 0)
start_dates_delta_vx[order(start_dates_delta_vx$Series_Complete_Pop_Pct_7, decreasing = TRUE), ]

# Mean Rt Estimates  ***********************************************************************

## Rt 
mean_rt_full_period <- ddply(all_rt, .(state, Variant_Category), summarize, Mean_Rt = round(mean(rt, na.rm = TRUE), digits = 2))
mean_rt_full_period_cast <- reshape2::dcast(mean_rt_full_period, state ~ Variant_Category, value.var = "Mean_Rt")
mean_rt_full_period_cast[order(mean_rt_full_period_cast$Other, decreasing = TRUE), ]
mean_rt_full_period_delta <- ddply(all_rt, .(Variant_Category), summarize, Mean_Rt = round(mean(rt, na.rm = TRUE), digits = 2))
mean_rt_full_period_delta %>% dplyr::filter(Variant_Category == "Delta")
mean_rt_full_period_cast[order(mean_rt_full_period_cast$Delta, decreasing = TRUE), ]

## Alpha mean Rt pre- and post-Delta
delta_dates <- start_end_dates_rbind %>% dplyr::filter(var == "Delta")
alpha_rt_periods <- all_rt %>% dplyr::filter(Variant_Category == "Alpha")
alpha_rt_periods$Period <- rep(NA, nrow(alpha_rt_periods))

for (state in unique(delta_dates$state)) {
  delta_start_date <- delta_dates[which(delta_dates$state == state), "min_date_var"]
  alpha_rt_periods[which((alpha_rt_periods$state == state) & (alpha_rt_periods$Date < delta_start_date)), "Period"] <- "Pre_Delta"
  alpha_rt_periods[which((alpha_rt_periods$state == state) & (alpha_rt_periods$Date >= delta_start_date)), "Period"] <- "Post_Delta"
}

alpha_pre_delt <- alpha_rt_periods %>% dplyr::filter(Period == "Pre_Delta")
alpha_pre_delt_mean_rt <- ddply(alpha_pre_delt, .(state), summarize, Mean_Rt = round(mean(rt, na.rm = TRUE), digits = 2))
alpha_pre_delt_mean_rt[order(alpha_pre_delt_mean_rt$Mean_Rt, decreasing = TRUE), ]
round(mean(alpha_pre_delt$rt, na.rm = TRUE), digits = 2)

alpha_post_delt <- alpha_rt_periods %>% dplyr::filter(Period == "Post_Delta")
alpha_post_delt_mean_rt <- ddply(alpha_post_delt, .(state), summarize, Mean_Rt = round(mean(rt, na.rm = TRUE), digits = 2))
alpha_post_delt_mean_rt[order(alpha_post_delt_mean_rt$Mean_Rt, decreasing = TRUE), ]
round(mean(alpha_post_delt$rt, na.rm = TRUE), digits = 2)

# Mean Rt Ratios  ***********************************************************************

mean_Rt_ratio <- ddply(delta_comp_all_limit, .(state, Comparison), summarize, Mean_Rt_Ratio = round(mean(Rt_Ratio, na.rm = TRUE), digits = 2))
mean_Rt_ratio_cast <- reshape2::dcast(mean_Rt_ratio, state ~ Comparison, value.var = "Mean_Rt_Ratio")
mean_Rt_ratio_cast[order(mean_Rt_ratio_cast$`Delta:Alpha`, decreasing = TRUE), ]

# Multiplicative Rt Increase  ***********************************************************************

mean_fold_inc <- ddply(model_df, .(lineage), summarize, Mean_Trans = round(mean(transmissibility), digits = 2))
mean_fold_inc_state <- ddply(model_df, .(state, lineage), summarize, Mean_Trans = round(mean(transmissibility), digits = 2))
mean_fold_inc_state_cast <- reshape2::dcast(mean_fold_inc_state, state ~ lineage, value.var = "Mean_Trans")
mean_fold_inc_state_cast[order(mean_fold_inc_state_cast$Delta, decreasing = TRUE), ]

