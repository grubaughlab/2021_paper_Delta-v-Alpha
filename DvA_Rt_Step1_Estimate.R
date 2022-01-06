#############################################################################################################################
# ABOUT ----------------------------------------------------------------------------------------------------
#############################################################################################################################

# This code generates the Rt estimates using a multi-step bootstrapping approach. The code is designed to be run on a cluster.

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
library(doParallel)

replicates <- 1000 
set.seed(1) # for reproducibility

# Load functions
path <- getwd()
file_name <- "/DvA_functions.R"
full_path <- paste0(path, file_name)
source(full_path, echo=TRUE)

# Load GISAID metadata file
metadata_date <- "2021-08-13"
path_data <- paste0(path, "data/")
load(paste0(path_data, "next_meta.Rda"))

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
next_meta_sub$Variant_Category <- as.character(next_meta_sub$Variant_Category)
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
# ESTIMATED INFECTIONS ----------------------------------------------
#############################################################################################################################

# https://pkg.covidestim.org/reference/summary.covidestim_result.html
# read in and format Covidestim data
file_name <- "covidestim_state_2021-11-22.csv"
infect_import <- read.csv(paste0(path_data, file_name))

infect = infect_import %>%
  dplyr::select(state, date, infections) %>% # median estimated infections/day, by date of infection
  dplyr::filter(state %in% states, date >= "2021-01-01", date <= "2021-07-31") %>% # matching dates for var frequencies (removing dates w/ NAs at start & end caused by 7 day rolling average)
  dplyr::rename(Date = date) %>%
  dplyr::mutate(Date = as.Date(Date))

infect_comb <- left_join(infect_old, infect_new, by = c("state" = 'state', "Date" = "Date"))
infect_comb$match <- infect_comb$infections.x == infect_comb$infections.y

#############################################################################################################################
#  DAILY VARIANT PROPORTIONS ----------------------------------------------
#############################################################################################################################

library(plyr)

# make sure all dates represented in variant frequency data 
date_range <- data.frame(rep(seq(as.Date("2021-01-01"), as.Date("2021-07-31"), by = "day"), length(states)))
colnames(date_range) <- "Date"
date_range$state <- rep(states, each = nrow(date_range) / 6)

next_states_melt_all <- next_states_melt %>% dplyr::filter(Date >= "2021-01-01", Date <= "2021-07-31")
next_states_melt_all_sum <- ddply(next_states_melt_all, .(state, Date, Variant_Category), summarize, Number_Genomes_Daily = sum(value))
next_states_melt_all_sum_dcast <- reshape2::dcast(next_states_melt_all_sum, state + Date ~ Variant_Category)

dat_all <- left_join(date_range, next_states_melt_all_sum_dcast, by = c("state" = "state", "Date" = "Date"))
dat_all[is.na(dat_all)] <- 0

# remove single Vermont outlier after dies out that heavily skews Rt
dat_all[which((dat_all$state == "Vermont") & (dat_all$Date == "2021-07-21")), "Other"] <- 0

# sample from multinomial distribution
state_samples <- list() 
for (j in states) { # for each state
  
  dat <- dat_all %>% dplyr::filter(state == j)
  dat <- dat[order(dat$Date, decreasing = FALSE), ]
  infect_state <- infect %>% dplyr::filter(state == j)

  var_probs <- ddply(dat, .(state, Date), summarize, 
                     Alpha = Alpha,
                     Delta = Delta,
                     Other = Other,
                     Num_Genomes_Daily = Alpha + Delta + Other,
                     Prop_Alpha = Alpha / Num_Genomes_Daily,
                     Prop_Delta = Delta / Num_Genomes_Daily,
                     Prop_Other = Other / Num_Genomes_Daily)
  
  var_probs[is.na(var_probs)] <- 0 
  
  s2_sampled_counts_all_days <- list()
  for (each_day in as.Date(unique(dat$Date))) {  # for each day
    
    # step 1
    s1_pop_day <- var_probs[which(var_probs$Date == each_day), "Num_Genomes_Daily"]
    s2_pop_day <- infect_state[which(infect_state$Date == each_day), "infections"]
    
    s1_var_probs_day <- var_probs[which(var_probs$Date == each_day), c("Prop_Alpha", "Prop_Delta", "Prop_Other")]
    s1_var_probs_day <- as.vector(s1_var_probs_day)
    
    if (s1_pop_day != 0) {
      s1_sampled_counts <- rmultinom(n = replicates, size = s1_pop_day, prob = s1_var_probs_day)
      s1_sampled_probs <- s1_sampled_counts / s1_pop_day
    }

    # step 2
    s2_sampled_counts <- list()
    each_day <- as.character(as.Date(each_day))
    
    for (each_rep in 1:replicates) { # for each replicate
      
      if (s1_pop_day != 0) { # if at least 1 genome on day i
        
        s2_var_probs_day <- s1_sampled_probs[, each_rep]
        s2_sampled_counts_single <- rmultinom(n = 1, size = s2_pop_day, prob = s2_var_probs_day)
        
        if (each_day == min(dat$Date)) {
          colnames(s2_sampled_counts_single) <- each_day
          s2_sampled_counts_all_days[[each_rep]] <- s2_sampled_counts_single
        } else {
          colnames(s2_sampled_counts_single) <- each_day
          s2_sampled_counts_all_days[[each_rep]] <- cbind.data.frame(s2_sampled_counts_all_days[[each_rep]], s2_sampled_counts_single)
        }
      } else { # if 0 genomes on day i
        
        if (each_day == min(dat$Date)) {
          s2_sampled_counts_single <- matrix(0, nrow = 3, ncol = 1)
          colnames(s2_sampled_counts_single) <- as.character(each_day)
          s2_sampled_counts_all_days[[each_rep]] <- s2_sampled_counts_single
        } else {
          s2_sampled_counts_single <- matrix(0, nrow = 3, ncol = 1)
          colnames(s2_sampled_counts_single) <- as.character(each_day)
          s2_sampled_counts_all_days[[each_rep]] <- cbind.data.frame(s2_sampled_counts_all_days[[each_rep]], s2_sampled_counts_single)
        }
      }
    }
  }
  state_samples[[j]] <- s2_sampled_counts_all_days
}

#############################################################################################################################
# EFFECTIVE REPRODUCTION NUMBER Rt ----------------------------------------------
#############################################################################################################################

# 7 day rolling average on estimated infections
state_samples_7 <- list()
for (j in names(state_samples)) {
  
  state_samples_1 <- state_samples[[j]]
  
  rep_list_7 <- list()
  for (each_rep in 1:length(state_samples_1)) {
    
    state_samples_1_rep <- state_samples_1[[each_rep]]
    state_samples_1_rep[1, ] <- zoo::rollmean(as.vector(data.matrix(state_samples_1_rep[1, ])), k = 7, fill = NA)
    state_samples_1_rep[2, ] <- zoo::rollmean(as.vector(data.matrix(state_samples_1_rep[2, ])), k = 7, fill = NA)
    state_samples_1_rep[3, ] <- zoo::rollmean(as.vector(data.matrix(state_samples_1_rep[3, ])), k = 7, fill = NA)
    state_samples_1_rep <- state_samples_1_rep[, 4:(ncol(state_samples_1_rep) - 3)] # get rid of NA cols induced by rolling avg
    rep_list_7[[each_rep]] <- state_samples_1_rep
  }
  state_samples_7[[j]] <- rep_list_7
}

save(state_samples_7, file = paste0(path_data, "state_samples_7", "_", replicates, ".rda"))  # need to have same folder name on cluster

# run Rt in parallel

## generate list of keys
iter_keys <- list()
ix <- 1
vars <- c("Prop_Alpha", "Prop_Delta", "Prop_Other")
reps <- 1:replicates
for (state in states) {
    for (r in reps) {
      for (variant in vars) {
      k <- list(state, r, variant)
      iter_keys[[ix]] <- k
      ix <- ix + 1
    }
  }
}

## function to return appropriate row matching to each key for input into EpiEstim
compute_rt <- function(k) {
  
  library(matrixStats)
  library(EpiEstim)
  
  state <- k[[1]]
  r <- k[[2]]
  variant <- k[[3]]
  
  # get data for this iteration key
  df <- state_samples_7[[state]][[r]][variant, ]
  
  # set so can get reproducible results
  MCMC_seed <- 1 
  overall_seed <- 2
  mcmc_control <- make_mcmc_control(seed = MCMC_seed)
  
  # calc Rt starting on day w/ at least 12 infections over 21 day windows
  output_rt <- rt_fun_replicate_parallel(df, window_size = 21, mcmc_control = mcmc_control, seed = overall_seed)
  
  list(list(k=k, result=output_rt))
} 

# setup parallel configuration
num_cores <-detectCores()
num_cores <- num_cores - 2
registerDoParallel(num_cores)

res <- foreach(k=iter_keys, .combine=c) %dopar% {
  print(k)
  compute_rt(k)
}

save(res, file = paste0(path_data, "res", "_", replicates, ".rda")) # need to have same folder name on cluster




