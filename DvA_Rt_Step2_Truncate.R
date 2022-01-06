#############################################################################################################################
# ABOUT ----------------------------------------------------------------------------------------------------
#############################################################################################################################

# This code checks the bootstapped Rt estimates, truncating variant-specific estimates on the date variant infections died out.

#############################################################################################################################
# LOAD LIBRARIES & FILES ----------------------------------------------------------------------------------------------------
#############################################################################################################################

library(RColorBrewer)
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)

# number of replicates in files
replicates <- 1000

# load functions
path <- getwd()
file_name <- "/DvA_functions.R"
full_path <- paste0(path, file_name)
source(full_path, echo=TRUE)

# path from where to grab data files
path_data <- paste0(path, "data/")

# specify selected states
states <- c("Connecticut", "Massachusetts", "Maine", "New Hampshire", "Rhode Island", "Vermont")

# select colors for plotting
all_dark2_colors <- brewer.pal(8, "Dark2")
customPalette <- c(all_dark2_colors[1], all_dark2_colors[3], all_dark2_colors[7]) # Alpha = green; Delta = purple; Other = brown

#############################################################################################################################
# DETERMINE TRUNCATION DATES ------------------------------------------------------------------------------------------------
#############################################################################################################################

# manually truncate Rt estimates plots where start to erroneously increase again despite all zeroes 
## (i.e. variant category died out due to mean prior of 5)
## only affects Alpha and Other

# load estimated variant-specific infections
load(paste0(path_data, "state_samples_7_", replicates, ".rda")) # this is output from Step 1

alpha_est_I <- list()
delta_est_I <- list()
other_est_I <- list()
dat <- state_samples_7
for (each_state in names(dat)) {
  
  dat_state <- dat[[each_state]]
  
  for (each_rep in 1:length(dat_state)) {
    
    dat_state_rep <- dat_state[[each_rep]]
    
    if (each_rep == 1) {
      alpha_state_all_reps <- data.frame(t(dat_state_rep[1, ]))
      delta_state_all_reps <- data.frame(t(dat_state_rep[2, ]))
      other_state_all_reps <- data.frame(t(dat_state_rep[3, ]))
    } else {
      alpha_state_all_reps <- cbind.data.frame(alpha_state_all_reps, data.frame(t(dat_state_rep[1, ])))
      delta_state_all_reps <- cbind.data.frame(delta_state_all_reps, data.frame(t(dat_state_rep[2, ])))
      other_state_all_reps <- cbind.data.frame(other_state_all_reps, data.frame(t(dat_state_rep[3, ])))
    }
    
  }
  
  alpha_state_all_reps$Date <- as.Date(rownames(alpha_state_all_reps))
  delta_state_all_reps$Date <- as.Date(rownames(delta_state_all_reps))
  other_state_all_reps$Date <- as.Date(rownames(other_state_all_reps))
  
  alpha_state_all_reps$state <- each_state
  delta_state_all_reps$state <- each_state
  other_state_all_reps$state <- each_state
  
  alpha_est_I[[each_state]] <- alpha_state_all_reps
  delta_est_I[[each_state]] <- delta_state_all_reps
  other_est_I[[each_state]] <- other_state_all_reps
}

alpha_est_I_all <- do.call(rbind.data.frame, alpha_est_I)
delta_est_I_all <- do.call(rbind.data.frame, delta_est_I)
other_est_I_all <- do.call(rbind.data.frame, other_est_I)

# fix naming
colnames(alpha_est_I_all) <- c(rep(paste0("Est_I_7_", seq(1:replicates))),"Date", "state")
colnames(delta_est_I_all) <- c(rep(paste0("Est_I_7_", seq(1:replicates))),"Date", "state")
colnames(other_est_I_all) <- c(rep(paste0("Est_I_7_", seq(1:replicates))),"Date", "state")

# format for plotting
selected_columns <- rep(paste0("Est_I_7_", seq(1:replicates)))
alpha_rt_melt <- reshape2::melt(alpha_est_I_all, id = c("Date", "state"), measure.vars = selected_columns)
delta_rt_melt <- reshape2::melt(delta_est_I_all, id = c("Date", "state"), measure.vars = selected_columns)
other_rt_melt <- reshape2::melt(other_est_I_all, id = c("Date", "state"), measure.vars = selected_columns)

# check die-outs (affects Rt estimates)
alpha_out_mean <- ddply(alpha_rt_melt, .(state, Date), summarize, Mean_Inf = mean(value, na.rm=TRUE))
alpha_out_mean <- alpha_out_mean[order(alpha_out_mean$Date, decreasing = TRUE), ]

# ## die-outs
# View(alpha_out_mean %>% filter(state == "Connecticut"))
# 2021-07-24
# View(alpha_out_mean %>% filter(state == "Massachusetts"))
# NA
# View(alpha_out_mean %>% filter(state == "Maine"))
# 2021-07-13
# View(alpha_out_mean %>% filter(state == "New Hampshire"))
# 2021-07-23
# View(alpha_out_mean %>% filter(state == "Rhode Island"))
# 2021-07-17
# View(alpha_out_mean %>% filter(state == "Vermont"))
# 2021-07-11

other_out_mean <- ddply(other_rt_melt, .(state, Date), summarize, Mean_Inf = mean(value, na.rm=TRUE))
other_out_mean <- other_out_mean[order(other_out_mean$Date, decreasing = TRUE), ]
# # ## die-outs
# View(other_out_mean %>% filter(state == "Connecticut"))
# # NA
# View(other_out_mean %>% filter(state == "Massachusetts"))
# # NA
# View(other_out_mean %>% filter(state == "Maine"))
# # NA
# View(other_out_mean %>% filter(state == "New Hampshire"))
# # 2021-07-20
# View(other_out_mean %>% filter(state == "Rhode Island"))
# # 2021-07-25
# View(other_out_mean %>% filter(state == "Vermont"))
# # 2021-06-25

# save manually selected truncation dates
alpha_ct <- "2021-07-24"
alpha_me <- "2021-07-13"
alpha_nh <- "2021-07-23"
alpha_ri <- "2021-07-17"
alpha_vt <- "2021-07-11"

other_nh <- "2021-07-20"
other_ri <- "2021-07-25"
other_vt <- "2021-06-25"

trunc_dates <- rbind.data.frame(alpha_ct, alpha_me, alpha_nh, alpha_ri, alpha_vt, other_nh, other_ri, other_vt)
colnames(trunc_dates) <- "Trunc_Date"
trunc_dates$state <- c("Connecticut", "Maine", "New Hampshire", "Rhode Island", "Vermont", "New Hampshire", "Rhode Island", "Vermont")
trunc_dates$Variant_Category <- c(rep("Alpha", 5), rep("Other", 3))
save(trunc_dates, file = paste0(path_data, "trunc_dates", "_", replicates, ".rda")) 

#############################################################################################################################
# PROCESS RAW RT DATA -------------------------------------------------------------------------------------------------------
#############################################################################################################################
load(paste0(path_data, "res_", replicates, ".rda")) # this is output from step 1

var_acc_list <- list()

for (entry in res) {
  k <- entry$k
  state <- k[[1]]
  r <- k[[2]]
  variant <- k[[3]]
  
  result <- entry$result
  colnames(result) <- c("Date", paste0("rt_", r))
  
  states_list <- var_acc_list[[variant]]
  if (is.null(states_list)) {
    states_list <- list()
  }
  
  replicates_df <- states_list[[state]]
  if (is.null(replicates_df)) {
    # nothing to do
    replicates_df <- result
  } else {
    # add new replicate as another column
    replicates_df <- left_join(replicates_df, result, by = c("Date" = "Date"))
  }
  # store it back
  states_list[[state]] <- replicates_df
  var_acc_list[[variant]] <- states_list
}  

# rbind and save truncated results together
alpha_rt <- add_state_rbind(var_acc_list$Prop_Alpha, states)
delta_rt <- add_state_rbind(var_acc_list$Prop_Delta, states)
other_rt <- add_state_rbind(var_acc_list$Prop_Other, states)

# add variant category
alpha_rt$Variant_Category <- "Alpha"
delta_rt$Variant_Category <- "Delta"
other_rt$Variant_Category <- "Other"

# take mean & 95% CI across states 
alpha_rt_mean_ci <- mean_CI_replicates(alpha_rt)
delta_rt_mean_ci <- mean_CI_replicates(delta_rt)
other_rt_mean_ci <- mean_CI_replicates(other_rt)

# combine into 1
all_rt_mean_ci <- rbind.data.frame(alpha_rt_mean_ci, delta_rt_mean_ci, other_rt_mean_ci)

#############################################################################################################################
# TRUNCATE MEAN RT ESTIMATES ------------------------------------------------------------------------------------------------
#############################################################################################################################

state_var_trunc_all <- NULL
dat <- all_rt_mean_ci
for (each_var in unique(trunc_dates$Variant_Category)) {
  
  trunc_dates_var <- trunc_dates %>% dplyr::filter(Variant_Category == each_var)
  
  state_trunc_all <- NULL
  for (each_state in states) {
    
    if (!(each_state%in% unique(trunc_dates_var$state))) {
      state_trunc_save <- dat %>% 
        dplyr::filter(state == each_state, Variant_Category == each_var)
    } else {
      trunc_date <- trunc_dates_var %>% 
        dplyr::filter(state == each_state) %>%
        dplyr::select("Trunc_Date")
      
      state_trunc_save <- dat %>% 
        dplyr::filter(state == each_state, Variant_Category == each_var, Date < as.Date(unlist(trunc_date)))
    }
    state_trunc_save$Variant_Category <- each_var
    state_trunc_all <- rbind.data.frame(state_trunc_all, state_trunc_save)
  }
  state_var_trunc_all <- rbind.data.frame(state_var_trunc_all, state_trunc_all)
}

delta_trunc_save <- dat %>% dplyr::filter(Variant_Category == "Delta")
all_rt_trunc_all <- rbind.data.frame(state_var_trunc_all, delta_trunc_save)
save(all_rt_trunc_all, file = paste0(path_data, "all_rt_trunc_all", "_", replicates, ".rda"))

#############################################################################################################################
# TRUNCATE RAW RT ESTIMATES ------------------------------------------------------------------------------------------------
#############################################################################################################################
# for ratio calculations 

## bind together
all_rt_raw <- rbind.data.frame(alpha_rt, delta_rt, other_rt)

state_var_trunc_all <- NULL
dat <- all_rt_raw
for (each_var in unique(trunc_dates$Variant_Category)) {
  
  trunc_dates_var <- trunc_dates %>% dplyr::filter(Variant_Category == each_var)
  
  state_trunc_all <- NULL
  for (each_state in states) {
    
    if (!(each_state%in% unique(trunc_dates_var$state))) {
      state_trunc_save <- dat %>% 
        dplyr::filter(state == each_state, Variant_Category == each_var)
    } else {
      trunc_date <- trunc_dates_var %>% 
        dplyr::filter(state == each_state) %>%
        dplyr::select("Trunc_Date")
      
      state_trunc_save <- dat %>% 
        dplyr::filter(state == each_state, Variant_Category == each_var, Date < as.Date(unlist(trunc_date)))
    }
    state_trunc_save$Variant_Category <- each_var
    state_trunc_all <- rbind.data.frame(state_trunc_all, state_trunc_save)
  }
  state_var_trunc_all <- rbind.data.frame(state_var_trunc_all, state_trunc_all)
}

delta_trunc_save <- dat %>% dplyr::filter(Variant_Category == "Delta")
all_rt_trunc_all_raw <- rbind.data.frame(state_var_trunc_all, delta_trunc_save)
save(all_rt_trunc_all_raw, file = paste0(path, "all_rt_trunc_all_raw", "_", replicates, ".rda"))

# # alpha 
# est_I_plot_alpha <- ggplot(alpha_rt_melt, aes(Date, value, group = variable)) +
#   geom_line(aes(color = variable)) +
#   # scale_y_continuous(breaks = seq(0, max_plot_mean, by = 1), limits = c(0, max_plot_mean)) +
#   facet_wrap(~state, ncol = 1, scales = "free_y") +
#   theme_bw() +
#   ylab("Number of Estimated Alpha Infections") +
#   theme(legend.position = "none") +
#   scale_x_date(date_breaks = "1 month")
# est_I_plot_alpha

# est_I_plot_alpha_end <- ggplot(alpha_rt_melt, aes(Date, value, group = variable)) +
#   geom_line(aes(color = variable)) +
#   facet_wrap(~state, ncol = 1) +
#   theme_bw() +
#   ylab("Number of Estimated Alpha Infections") +
#   theme(legend.position = "none") +
#   scale_x_date(date_breaks = "1 month", limits = c(as.Date("2021-05-01"), NA)) +
#   ylim(0, 30)
# est_I_plot_alpha_end
# 
# # delta
# est_I_plot_delta <- ggplot(delta_rt_melt, aes(Date, value, group = variable)) +
#   geom_line(aes(color = variable)) +
#   # scale_y_continuous(breaks = seq(0, max_plot_mean, by = 1), limits = c(0, max_plot_mean)) +
#   facet_wrap(~state, ncol = 1, scales = "free_y") +
#   theme_bw() +
#   ylab("Number of Estimated Delta Infections") +
#   theme(legend.position = "none") +
#   scale_x_date(date_breaks = "1 month")
# est_I_plot_delta
# 
# # other
# est_I_plot_other <- ggplot(other_rt_melt, aes(Date, value, group = variable)) +
#   geom_line(aes(color = variable)) +
#   # scale_y_continuous(breaks = seq(0, max_plot_mean, by = 1), limits = c(0, max_plot_mean)) +
#   facet_wrap(~state, ncol = 1, scales = "free_y") +
#   theme_bw() +
#   ylab("Number of Estimated Other Infections") +
#   theme(legend.position = "none") +
#   scale_x_date(date_breaks = "1 month")
# est_I_plot_other
# 
# est_I_plot_other_end <- ggplot(other_rt_melt, aes(Date, value, group = variable)) +
#   geom_line(aes(color = variable)) +
#   theme_bw() +
#   ylab("Number of Estimated Other Infections") +
#   theme(legend.position = "none") +
#   scale_x_date(date_breaks = "1 month", limits = c(as.Date("2021-05-01"), NA)) +
#   facet_wrap(~state, ncol = 1, scales = "free_y") +
#   xlim(c(as.Date("2021-06-01"), NA))  +
#   ylim(0, 30)
# est_I_plot_other_end
