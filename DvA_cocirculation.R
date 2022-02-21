#############################################################################################################################
# ABOUT ----------------------------------------------------------------------------------------------------
#############################################################################################################################

# This code calculates Alpha and Delta co-circulation periods in each state.
## Figure S3 and Table S5 in the supplementary text

#############################################################################################################################
# LOAD LIBRARIES & FILES ----------------------------------------------------------------------------------------------------
#############################################################################################################################

library(RColorBrewer)
library(ggplot2)
library(scales)
library(dplyr)
library(plyr)

# number of replicates in files
replicates <- 1000

# load functions
path <- getwd()
file_name <- "/DvA_functions.R"
full_path <- paste0(path, file_name)
source(full_path, echo=TRUE)

# Load GISAID metadata file 
metadata_date <- "2021-08-13"
gisaid_date <- paste("/", metadata_date, "/", sep="")
path_data <- paste0(path, "data/")
load(paste0(path_data, "/next_meta.Rda"))

# specify selected states
states <- c("Connecticut", "Massachusetts", "Maine", "New Hampshire", "Rhode Island", "Vermont")

# select colors for plotting
all_dark2_colors <- brewer.pal(8, "Dark2")
customPalette <- c(all_dark2_colors[1], all_dark2_colors[3], all_dark2_colors[7]) # Alpha = green; Delta = purple; Other = brown

# path for where to save figures
path_figures <- paste0(path, "figures/")

#############################################################################################################################
# Delta and Alpha Co-Circulation - Infections  ----------------------------------------------
#############################################################################################################################

# load estimated variant-specific infections
load(paste0(path_data, "state_samples_7_", replicates, ".rda")) # output from Rt estimation step 1

alpha_est_I <- list()
delta_est_I <- list()
dat <- state_samples_7
for (each_state in names(dat)) {
  
  dat_state <- dat[[each_state]]
  
  for (each_rep in 1:length(dat_state)) {
    
    dat_state_rep <- dat_state[[each_rep]]
    
    if (each_rep == 1) {
      alpha_state_all_reps <- data.frame(t(dat_state_rep[1, ]))
      delta_state_all_reps <- data.frame(t(dat_state_rep[2, ]))
    } else {
      alpha_state_all_reps <- cbind.data.frame(alpha_state_all_reps, data.frame(t(dat_state_rep[1, ])))
      delta_state_all_reps <- cbind.data.frame(delta_state_all_reps, data.frame(t(dat_state_rep[2, ])))
    }
    
  }
  
  alpha_state_all_reps$Date <- as.Date(rownames(alpha_state_all_reps))
  delta_state_all_reps$Date <- as.Date(rownames(delta_state_all_reps))

  alpha_state_all_reps$state <- each_state
  delta_state_all_reps$state <- each_state

  alpha_est_I[[each_state]] <- alpha_state_all_reps
  delta_est_I[[each_state]] <- delta_state_all_reps
}

alpha_est_I_all <- do.call(rbind.data.frame, alpha_est_I)
delta_est_I_all <- do.call(rbind.data.frame, delta_est_I)

# fix naming
colnames(alpha_est_I_all) <- c(rep(paste0("Est_I_7_", seq(1:replicates))),"Date", "state")
colnames(delta_est_I_all) <- c(rep(paste0("Est_I_7_", seq(1:replicates))),"Date", "state")

# format for plotting
selected_columns <- rep(paste0("Est_I_7_", seq(1:replicates)))
alpha_I_melt <- reshape2::melt(alpha_est_I_all, id = c("Date", "state"), measure.vars = selected_columns)
delta_I_melt <- reshape2::melt(delta_est_I_all, id = c("Date", "state"), measure.vars = selected_columns)

# calculate alpha die-out date for each state, replicate
die_out_date_list <- list()
dat <- alpha_I_melt
for (each_state in unique(dat$state)) {
  
  dat_state <- dat %>% dplyr::filter(state == each_state)
  
  die_out_state <- NULL
  for (each_rep in unique(dat_state$variable)) {
    dat_state_rep <- dat_state %>% dplyr::filter(variable == each_rep)
    dat_state_rep <- dat_state_rep[order(dat_state_rep$Date, decreasing = TRUE), ]
    dat_state_rep[which(dat_state_rep$value == 0), "Zero"] <- 0
    dat_state_rep[which(dat_state_rep$value != 0), "Zero"] <- 1
    dat_state_rep$Zero_CumSum <- cumsum(dat_state_rep$Zero)
    die_out <- ifelse(nrow(dat_state_rep[which(dat_state_rep$Zero_CumSum == 0), ]) == 0, 
                      max(dat_state_rep$Date), 
                      min(dat_state_rep[which(dat_state_rep$Zero_CumSum == 0), "Date"]))
    die_out_state_single <- cbind.data.frame(Rep = each_rep, Die_Out_Date = as.Date(die_out, origin = "1970-01-01"))
    die_out_state <- rbind.data.frame(die_out_state, die_out_state_single)
  }
  die_out_date_list[[each_state]] <- die_out_state
}

# calculate delta start date for each state, replicate
start_date_list <- list()
dat <- delta_I_melt
for (each_state in unique(dat$state)) {
  
  dat_state <- dat %>% dplyr::filter(state == each_state)
  
  start_state <- NULL
  for (each_rep in unique(dat_state$variable)) {
    dat_state_rep <- dat_state %>% dplyr::filter(variable == each_rep)
    dat_state_rep <- dat_state_rep[order(dat_state_rep$Date, decreasing = FALSE), ]
    min_start_date <- min(dat_state_rep[which(dat_state_rep$value > 0), "Date"])
    start_state_single <- cbind.data.frame(Rep = each_rep, Start_Date = as.Date(min_start_date))
    start_state <- rbind.data.frame(start_state, start_state_single)
  }
  start_date_list[[each_state]] <- start_state
}

# for each state, replicate count number of co-circulation days
mean_cocirc_inf_list <- list()
for (each_state in names(die_out_date_list)) {

  die_state <- die_out_date_list[[each_state]]
  start_state <- start_date_list[[each_state]]
  
  start_die_state <- left_join(start_state, die_state, by = c("Rep" = "Rep"))
  start_die_state$Diff <- start_die_state$Die_Out_Date - start_die_state$Start_Date
  mean_start_inf <- mean(start_die_state$Start_Date)
  mean_die_inf <- mean(start_die_state$Die_Out_Date)
  mean_cocirc_inf <- mean(start_die_state$Diff)
  state_cocirc_inf <- cbind.data.frame(State = each_state, Mean_Start_Date = mean_start_inf, Mean_Die_Date = mean_die_inf, Mean_Cocirc_Pd = mean_cocirc_inf)
  mean_cocirc_inf_list[[each_state]] <- state_cocirc_inf
}

mean_cocirc_inf_all <- do.call(rbind.data.frame, mean_cocirc_inf_list)
mean_cocirc_inf_all$Mean_Cocirc_Pd <- round(mean_cocirc_inf_all$Mean_Cocirc_Pd, digits = 0)
rownames(mean_cocirc_inf_all) <- NULL
colnames(mean_cocirc_inf_all) <- c("State", "Mean Delta Start Date" , "Mean Alpha Die-Out Date", "Mean Co-Circulation Period")
mean_cocirc_inf_all <- mean_cocirc_inf_all[order(mean_cocirc_inf_all$State, decreasing = FALSE), ]

library(kableExtra)
cocirc_inf <- kableExtra::kable(mean_cocirc_inf_all) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, position = "left")
# screenshot and save as pdf - Table_S5

# plot mean alpha and delta infections across replicates
alpha_I_melt$Variant_Category <- "Alpha"
delta_I_melt$Variant_Category <- "Delta"
alpha_delta_I <- rbind.data.frame(alpha_I_melt, delta_I_melt)

# restrict plot by mean start and die-out dates for each state
dat <- alpha_delta_I
dat_cocirc <- mean_cocirc_inf_all
trunc_inf_list <- list()
for (each_state in unique(dat$state)) {
  
  dat_state <- dat %>% dplyr::filter(state == each_state)
  dat_cocirc_state <- dat_cocirc %>% 
    dplyr::filter(State == each_state) %>%
    dplyr::select("Mean Delta Start Date", "Mean Alpha Die-Out Date")
  earliest_date <- as.Date(unlist(dat_cocirc_state[1]), origin = "1970-01-01")
  latest_date <-  as.Date(unlist(dat_cocirc_state[2]), origin = "1970-01-01")
  dat_state_dates <- dat_state %>% dplyr::filter(Date >= earliest_date, Date <= latest_date)
  trunc_inf_list[[each_state]] <- dat_state_dates
}

trunc_inf_all <- do.call(rbind.data.frame, trunc_inf_list)

# add 2019 state populations
dat_pop <- read_excel("nst-est2019-01.xlsx")
dat_pop <- dat_pop[4:nrow(dat_pop), c(1, 13)]
colnames(dat_pop) <- c("state", "pop")
dat_pop$state <- gsub("[.]", "", dat_pop$state)
dat_pop <- dat_pop %>% dplyr::filter(state %in% states)

trunc_inf_all_pop <- left_join(trunc_inf_all, dat_pop, by  = c("state" = "state"))
trunc_inf_all_pop$Inf_per_100K <- (trunc_inf_all_pop$value / trunc_inf_all_pop$pop) * 100000
trunc_inf_all_pop_mean <- ddply(trunc_inf_all_pop, .(state, Variant_Category, Date), summarize, Mean_Inf_per_100K = mean(Inf_per_100K, na.rm = TRUE))

p <- ggplot(trunc_inf_all_pop_mean) + 
  geom_bar(aes(x = Date, y = Mean_Inf_per_100K, fill = Variant_Category, color = Variant_Category), stat = "identity") +
  scale_fill_manual(values = customPalette, name = "Variant Category") +
  scale_color_manual(values = customPalette, name = "Variant Category") +
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")) + 
  facet_wrap(~state, scales = "fixed") +
  theme(legend.position="bottom") +
  ylab("Mean Estimated Infections per 100K Population") 

save(trunc_inf_all_pop_mean, file = paste0(path_data, "cocirc_plot_data.rda"))

#############################################################################################################################
# Delta and Alpha Co-Circulation - Rt  ----------------------------------------------
#############################################################################################################################

load(paste0(path_data, "all_rt_trunc_all_", replicates, ".rda"))

# calculate last Alpha Rt estimates for each state
all_rt_trunc_all_val <- all_rt_trunc_all %>% dplyr::filter(is.na(mean) == FALSE)
max_date <- ddply(all_rt_trunc_all_val, .(state, Variant_Category), summarize, Max_Date = max(Date))
max_date_alpha <- max_date %>% dplyr::filter(Variant_Category == "Alpha")

# calculate 1st Delta Rt estimates for each state
min_date <- ddply(all_rt_trunc_all_val, .(state, Variant_Category), summarize, Min_Date = min(Date))
min_date_delta <- min_date %>% dplyr::filter(Variant_Category == "Delta")

# combine min Delta and max Alpha dates
combo_dates <- left_join(min_date_delta, max_date_alpha, by = c("state" = "state"))
combo_dates$Cocirc_Time <- combo_dates$Max_Date - combo_dates$Min_Date

# save co-circulation dates
cocirc_dates <- combo_dates[, c("state", "Min_Date", "Max_Date")]
colnames(cocirc_dates) <- c("state", "Delta_Start", "Alpha_End")
save(cocirc_dates, file = paste0(path_data, "cocirc_dates.rda"))