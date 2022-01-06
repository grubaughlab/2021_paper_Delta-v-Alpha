#############################################################################################################################
# ABOUT ----------------------------------------------------------------------------------------------------
#############################################################################################################################

# This code performs a sensitivity analysis of the emergence period used in calculating the logistic growth rates.  
## Figure S2 and Table S4 in the supplementary material

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
library(MASS)

# Load functions
path <- getwd()
file_name <- "/DvA_functions.R"
full_path <- paste0(path, file_name)
source(full_path, echo=TRUE)

# Load GISAID metadata file 
metadata_date <- "2021-08-13"
gisaid_date <- paste("/", metadata_date, "/", sep="")
path_data <- paste0(path, "data/")
load(paste0(path_data, "/next_meta.Rda"))

# set the number of days for the emergence period
emergence_pd_vec <- c(60, 90, 120)

# path for where to save figures
path_figures <- paste0(path, "figures/")

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
# LOGISTIC GROWTH DURING VARIANT-SPECIFIC EMERGENCE PERIODS ----------------------------------------------
#############################################################################################################################
coef_bar_emerge_list <- list()
all_coefs_list <- list()
growth_plots_emerge_list <- list()
var_categories <- c("Alpha", "Delta")

for (emergence_pd in emergence_pd_vec) {
  # calculate start and end dates of 90 day emergence periods for each variant category / state
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
  growth_plots_emerge_list[[as.character(emergence_pd)]] <- plotLineage_cat_Xdays_multiple(next_states_melt_emerge, num_days = emergence_pd)
  
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
  
  # plot binomial logistic regression coefficients bar chart w/ CIs
  all_coefs <- do.call(rbind.data.frame, state_coefs)
  all_coefs <- all_coefs[, c("lineage", "coefs", "lowci", "upci")]
  all_coefs$state <- rep(states, each = length(unique(all_coefs$lineage)))
  all_coefs$log_odds <- all_coefs$coefs
  all_coefs$lineage <- factor(all_coefs$lineage, levels = c("Delta", "Alpha"))
  all_coefs_list[[as.character(emergence_pd)]] <- all_coefs
  coef_bar_emerge_list[[as.character(emergence_pd)]] <- plot_coefs_bar_sep_rotate(all_coefs)
}

#############################################################################################################################
# Figure S2: Slope for Alpha vs Delta over Differing Emergence Periods ----------------------------------------------
#############################################################################################################################
library(gridExtra)
library(ggpubr)

all_coefs_combined <- NULL
for (each_period in names(all_coefs_list)) {
  
  all_coefs_single <- all_coefs_list[[each_period]]
  all_coefs_single$period <- rep(each_period, nrow(all_coefs_single))
  all_coefs_combined <- rbind.data.frame(all_coefs_combined, all_coefs_single)
}

all_coefs_combined$lineage <- factor(all_coefs_combined$lineage, 
                                    levels = c("Alpha", "Delta"))
all_coefs_combined$period <- factor(all_coefs_combined$period, 
                                    levels = c("120", "90", "60"))
all_coefs_combined$state <- factor(all_coefs_combined$state, 
                                    levels = rev(states))

customPalette_2 <- brewer.pal(8, "RdBu") 

all_coefs_combined_alpha <- all_coefs_combined %>% dplyr::filter(all_coefs_combined$lineage == "Alpha")
all_coefs_combined_delta <- all_coefs_combined %>% dplyr::filter(all_coefs_combined$lineage == "Delta")

all_coefs_combined_alpha$state <- factor(all_coefs_combined_alpha$state, levels = states[order(states, decreasing = TRUE)])
all_coefs_combined_delta$state <- factor(all_coefs_combined_delta$state, levels = states[order(states, decreasing = TRUE)])

p_alpha <- plot_coefs_bar_sep_rotate_multiple(all_coefs_combined_alpha) + scale_fill_manual(values = customPalette_2)
p_delta <- plot_coefs_bar_sep_rotate_multiple(all_coefs_combined_delta) + scale_fill_manual(values = customPalette_2)

mylegend <- g_legend(p_alpha  + theme(legend.position='bottom'))
p <- grid.arrange(ggarrange(p_alpha + theme(legend.position="none"),
                            p_delta + theme(legend.position="none"),
                            ncol = 2, labels = c("A", "B")),
                  mylegend,
                  nrow = 2,
                  heights = c(9.5, 0.5))

ggsave(paste(path_figures, "Fig_S2.pdf", sep=""), p, height = 10, width =10)

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
# add relative logistic regression coefficients: Delta / Alpha
all_coefs_cast <- reshape2::dcast(all_coefs_combined, state + period ~ lineage, value.var = "log_odds")
all_coefs_cast$delta_div_alpha_perc <- ((all_coefs_cast$Delta - all_coefs_cast$Alpha)/all_coefs_cast$Alpha) * 100
all_coefs_cast$delta_div_alpha <- (all_coefs_cast$Delta/all_coefs_cast$Alpha)

all_coefs_cast$delta_div_alpha_perc <- round(all_coefs_cast$delta_div_alpha_perc, digits = 0)
all_coefs_cast$delta_div_alpha <- round(all_coefs_cast$delta_div_alpha, digits = 2)

all_coefs_cast[order(all_coefs_cast$delta_div_alpha, decreasing = TRUE), ] 

all_coefs_cast_nh_me_60 <- all_coefs_cast %>% dplyr::filter(state %in% c("New Hampshire", "Maine"), period == 60)

all_coefs_cast_2 <- reshape2::dcast(all_coefs_cast, state ~ period, value.var = "delta_div_alpha")
all_coefs_cast_2 <- all_coefs_cast_2[, c("state", "60", "90", "120")]
colnames(all_coefs_cast_2) <- c("State", "60", "90", "120")
all_coefs_cast_2$State <- as.character(all_coefs_cast_2$State)
all_coefs_cast_2 <- all_coefs_cast_2[order(all_coefs_cast_2$State, decreasing = FALSE), ]

rownames(all_coefs_cast_2) <- NULL
library(kableExtra)
growth_times <- kableExtra::kable(all_coefs_cast_2, format = "html") %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left") %>% 
  add_header_above(c(" " = 1, "Emergence Period (Days)" = 3))

# save_kable(growth_times, paste0(path_figures, "Table_S4.html"))
# screenshot & save as pdf

all_coefs_cast_2[order(all_coefs_cast_2$`60`, decreasing = TRUE), ]
all_coefs_cast_2[order(all_coefs_cast_2$`120`, decreasing = TRUE), ]

