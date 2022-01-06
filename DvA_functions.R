
# process initial GISAID dataset
readNextMeta_states <- function(next_meta_path){
  SS <- read_tsv(next_meta_path)
  SS <- SS[,c("Virus name", "Collection date", "Pango lineage", "AA Substitutions","Clade", "Location", "Accession ID")]
  names(SS) <- c("seqName", "Date", "pango_lineage", "AAsubstitutions", "Clade", "Location", "Accession ID")
  SS$Date <- ymd(SS$Date)
  SS$seqName <- gsub("hCoV-19/", "", SS$seqName)
  SS$country <- sapply(SS$seqName, function(x) str_split(x, "/")[[1]][1])
  SS$state_1 <- substr(SS$seqName, 5,6)
  SS$state_2 <-  sapply(SS$Location, function(x) str_split(x, "/")[[1]][3])
  SS$epiweek <- epiweek(SS$Date)
  SS <- SS %>% mutate(month = month(Date)) %>%
    mutate(year = year(Date)) %>%
    mutate(week = week(Date)) %>%
    filter(!is.na(year) & !is.na(month)) 
  SS <- convertMonth(SS)
  SS$AAsubstitutions <- gsub("[()]", "", SS$AAsubstitutions)
  SS$Sgeno <- sapply(SS$AAsubstitutions, extractS2) # parse nextclade format to extract S genotype
  SS$Sgeno <- gsub("Spike_", "S:", SS$Sgeno)
  SS$epidate <- get_date(week = SS$week, year = SS$year) # gives 1st day of isoweeks
  SS
}

extractS2 <- function(mutlist){
  splitlist <- strsplit(mutlist, ",")[[1]]
  paste(splitlist[grep("Spike_", splitlist)], collapse=",")
}

convertMonth <- function(SS){
  SS$jmonth <- SS$month + (12*(SS$year - 2020))
  SS$jweek <- SS$week + (53*(SS$year - 2020))
  SS
}


# label each variant category w/ binary 0/1
encodeLineage_var <- function(SS){
  one_hot_encoding <- one_hot(as.data.table(as.factor(SS$Variant_Category)))
  names(one_hot_encoding) <- gsub("V1_", "", names(one_hot_encoding))
  cbind(SS, one_hot_encoding)
}

# calculate date of 1st variant category detection in each state
calc_date_1st_inf <- function(next_states_melt, var_categories, num_days) {
  
  dat <- next_states_melt
  
  state_list <- list()
  for (j in states) {
    
    df_state <- dat %>% dplyr::filter(state == j)
    
    df_all <- NULL
    for (var in var_categories) {
      df_state_var <- df_state[which(df_state$Variant_Category == var), ]
      df_state_var <- df_state_var[which(df_state_var$value == 1), ]
      min_date_var <- as.Date(min(df_state_var$Date))
      max_date_var <- as.Date(max(seq.Date(from = min_date_var, length.out = num_days, by = "day")))
      df_single <- cbind.data.frame(var, min_date_var, max_date_var)
      df_all <- rbind.data.frame(df_all, df_single)
    }
    state_list[[j]] <- df_all
  }
  state_list
}

# restrict data to X day emergence periods; add counter
transform_data_Xdays <- function(dat, start_end_dates, num_days) {
  
  state_list <- list()
  for (j in states) {
    df_state <- dat %>% dplyr::filter(state == j)
    dates_state <- start_end_dates[[j]]
    
    df_all <- NULL
    for (var in var_categories) {
      df_state_var <- df_state[which(df_state$Variant_Category == var), ]
      df_dates_var <- dates_state[which(dates_state$var == var), ]
      my_dates <- seq.Date(df_dates_var$min_date_var, length.out = num_days, by = "day")
      df_state_var_date <- df_state_var[which((df_state_var$Date >= min(my_dates)) & (df_state_var$Date  <= max(my_dates))), ]
      dates_ticker <- cbind.data.frame(my_dates, 1:num_days)
      colnames(dates_ticker) <- c("Date", "Counter")
      df_single <- left_join(df_state_var_date, dates_ticker, by = c("Date" = "Date"))
      df_all <- rbind.data.frame(df_all, df_single)
    }
    state_list[[j]] <- df_all
  }
  state_list
}

# plot logistic growth curves
plotLineage_cat_Xdays <- function(SS, num_days){
  alphapoint = 0.5
  alphaline = 0.8
  alphashade = 0.2
  ggplot(SS, aes(x = Counter, y = value, color = Variant_Category, fill = Variant_Category)) + 
    theme_bw()+ 
    scale_fill_manual(values = customPalette) +
    scale_color_manual(values = customPalette) +
    geom_line(method = "glm",method.args=list(family="binomial"), alpha = alphaline, fullrange=TRUE, stat="smooth") +
    stat_smooth(method = "glm",method.args=list(family="binomial"), alpha = alphashade, fullrange=TRUE, size = 0) +
    facet_wrap(~state, nrow = length(states), ncol = 1) +
    labs(color = "Variant Category", fill = "Variant Category") + 
    ylim(0, 1) +
    scale_x_continuous(breaks = seq(0, num_days, by = 15)) +
    xlab("Days Since 1st State Variant Sequence") +
    ylab("Probability of Belonging to Variant Category") +
    theme(legend.position = "bottom", legend.box = "horizontal")
}

# plot binomial logistic regression coefficients bar chart w/ CIs
plot_coefs_bar_sep_rotate <- function(df){
  ggplot(df, aes(x = lineage, y = log_odds, fill = lineage)) + 
    theme_bw()+ 
    geom_bar(stat='identity', position = "dodge") +
    geom_errorbar(aes(ymin=lowci, ymax=upci), width = 0.2) +
    scale_color_manual(values = customPalette[2:1]) +  
    scale_fill_manual(values = customPalette[2:1]) +
    ylab("Slope") + # Log Odds of Belonging to Each Variant Category
    xlab("Variant Category") + 
    facet_wrap(~state, nrow = 6, ncol = 1) +
    coord_flip()
}

# calculate daily proportions of each variant category (note: days w/ no sequences will not be represented)
calc_prop_day <- function(dat, var_categories, timeInt){
  
  state_keep <- list()
  for (lineage in var_categories) {
    l_int <- dat %>% 
      dplyr::group_by(.data[[timeInt]]) %>% 
      dplyr::summarise(Proportion = mean(.data[[lineage]], na.rm=T), K= sum(.data[[lineage]], na.rm=T), N = n())
    state_keep[[lineage]] <- l_int
  }
  state_keep
}

# check for long stretches of zero variant-specific infections - can affect Rt estimates 
check_outliers <- function(variant_df_vec, num_zero_days, interest_col) {
  
  prob_log_list <- list()
  
  for (each_state in states) {
    
    dat <- variant_df_vec %>% dplyr::filter(state == each_state)
    dat_I <- dat[, interest_col]
    dat_I[which(is.na(dat_I))] <- 0
    dat_I <- dat_I[(min(which(dat_I > 0))):length(dat_I)]
    if ((any(rle(dat_I)$lengths[rle(dat_I)$values==0] > num_zero_days)) == TRUE) {
      prob_log_list[[state]] <- state
    }
  }
  prob_log <- do.call(rbind.data.frame, prob_log_list)
  prob_log
}

# check for long stretches of zero variant-specific infections - can affect Rt estimates 
check_outliers_replicate <- function(variant_df_vec, num_zero_days, interest_col) {
  
  prob_log_list <- list()
  
  for (each_state in states) {
    
    dat <- variant_df_vec %>% dplyr::filter(state == each_state)
    
    length_zeroes <- NULL
    prob_log_reps <- NULL
    for (each_rep in 1:length(unique(dat$variable))) { # for each replicate column
      
      interest_col_rep <- paste0(interest_col, each_rep) 
      dat_I <- dat[which(dat$variable == interest_col_rep), ]
      dat_I <- dat_I[order(dat_I$Date, decreasing = FALSE), ]
      dat_I[which(is.na(dat_I))] <- 0
      dat_I <- dat_I[min(which(dat_I$value > 0)):nrow(dat_I), ] # gets rid of NA/0s at start
      
      # if we have any 0 values that appears in a row for more than the number of days
      if ((any(rle(dat_I$value)$lengths[rle(dat_I$value)$values==0] > num_zero_days)) == TRUE) {
        prob_log_reps <- c(prob_log_reps, each_rep)
        max_dieout_rep <- rle(dat_I$value)$lengths[rle(dat_I$value)$values==0] # which died out
        max_dieout_rep <- max_dieout_rep[max_dieout_rep > num_zero_days] # max of dieout days > num_zero_days
        max_dieout_rep <- max(max_dieout_rep)
        length_zeroes <- c(length_zeroes, max_dieout_rep)
      }
    }
    state_results <- rbind.data.frame(prob_log_reps, length_zeroes)
    if (length(state_results != 0)) {
      colnames(state_results) <- NULL
      rownames(state_results) <- c("Replicate_Number", "Max_Dieout_Length")
    }
    prob_log_list[[each_state]] <- state_results
  }
  prob_log_list
}

# rbinds list of states and adds name
add_state_rbind <- function(dat, states) {
  
  dat_all <- NULL
  for (state in states) {
    dat_state <- dat[[state]]
    dat_state$state <- rep(state, nrow(dat_state))
    dat_all <- rbind.data.frame(dat_all, dat_state)
  }
  dat_all
}

# plot vaccinations vs infections
plot_vx_inf <- function(dat_vx_inf_pop, start_end_dates_rbind) {
  
  dat <- dat_vx_inf_pop
  max_y <- max(dat$Inf_per_100K, na.rm = TRUE) 
  
  sec_axis_scale <-  max(dat$Inf_per_100K, na.rm = TRUE) / 100 # b/c need extra space to plot labels
  
  p <- ggplot(dat) + 
    geom_bar(mapping = aes(x = Date, y = Inf_per_100K), fill = "grey", stat = "identity", position = "stack") +
    geom_line(mapping = aes(x = Date, y = Series_Complete_Pop_Pct_7*sec_axis_scale)) + 
    scale_color_manual(values = customPalette) +  
    scale_fill_manual(values = customPalette) +    
    scale_y_continuous(name = "Estimated Infections per 100K Population",
                       sec.axis = sec_axis(~./sec_axis_scale, name = "Percent of Population Fully Vaccinated",
                                           labels = function(b) { paste0(round(b, 0), "%")})) +
    geom_rect(data = start_end_dates_rbind, aes(xmin =  min_date_var, xmax = max_date_var, ymin = -Inf, ymax = max_y, fill = var), alpha = 0.15) +
    facet_wrap(~state, nrow = 6, ncol = 1) +
    labs(fill = "Variant Category") +
    xlab("Date")  +
    theme_bw() +
    theme(legend.position = "bottom", 
          axis.title.y = element_text(color = "black"),
          axis.title.y.right = element_text(color = "black")) 
  p
}


# plot proportions over time
plotProp_line <- function(SS, var_categories, timeInt){
  
  df <- NULL
  for (lineage in var_categories) {
    
    l_int <- SS %>% 
      group_by(state, epidate) %>% # for each state and week, perform the below
      dplyr::summarise(Proportion = mean(.data[[lineage]], na.rm=T), K= sum(.data[[lineage]], na.rm=T), N = n())
    
    l_int$LowerCI <- BinomCI(l_int$K, l_int$N, method = "jeffreys")[,2] # calculate 95% Jeffreys' interval
    l_int$UpperCI <- BinomCI(l_int$K, l_int$N, method = "jeffreys")[,3]
    l_int$Variant_Category <- rep(lineage, nrow(l_int))
    df <- rbind.data.frame(df, l_int)
  }
  
  df$Variant_Category <- factor(df$Variant_Category, levels = c("Alpha", "Delta","Other")) # converts to factor format for plotting 
  
  pp <- ggplot(df, aes(x = epidate, y = Proportion, color = Variant_Category, fill = Variant_Category)) + 
    geom_line() +
    geom_ribbon(aes(ymin=LowerCI, ymax=UpperCI), linetype = 0, alpha = 0.2) +
    scale_fill_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) +
    scale_color_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) +
    # geom_smooth(method = "loess", se = TRUE) +
    theme_bw() + 
    xlab("Date") +
    labs(color = "Variant Category", fill = "Variant Category") + 
    theme(legend.position="bottom") +
    facet_wrap(~state, ncol = 1) +
    scale_y_continuous(limits = c(0, 1)) 
  pp
}

# for grabbing legend for use in a combination plot
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mean_CI_replicates <- function(dat) {
  
  dat_row_all <- NULL
  for (each_row in 1:nrow(dat)) {
    
    dat_row <- dat[each_row, ]
    
    if(any(!is.na(dat_row[, 2:(1+replicates)]))) {
      dat_row$mean <- rowMeans(dat_row[, 2:(1+replicates)], na.rm = TRUE)
      dat_row$loci <- quantile(dat_row[, 2:(1+replicates)], prob = 0.025, na.rm = TRUE)
      dat_row$hici <- quantile(dat_row[, 2:(1+replicates)], prob = 0.975, na.rm = TRUE) 
    } else {
      dat_row$mean <- NA
      dat_row$loci <- NA
      dat_row$hici <- NA 
    }
    dat_row <- dat_row[, c("state", "Date", "Variant_Category", "mean", "loci", "hici")]
    dat_row_all <- rbind.data.frame(dat_row_all, dat_row)
  }
  dat_row_all
}

mean_CI_replicates_ratio <- function(dat) {
  
  dat_row_all <- NULL
  for (each_row in 1:nrow(dat)) {
    
    dat_row <- dat[each_row, ]
    
    if(any(!is.na(dat_row[, 3:(1+replicates)]))) {
      dat_row$mean <- rowMeans(dat_row[, 3:(2+replicates)], na.rm = TRUE)
      dat_row$loci <- quantile(dat_row[, 3:(2+replicates)], prob = 0.025, na.rm = TRUE)
      dat_row$hici <- quantile(dat_row[, 3:(2+replicates)], prob = 0.975, na.rm = TRUE) 
    } else {
      dat_row$mean <- NA
      dat_row$loci <- NA
      dat_row$hici <- NA 
    }
    dat_row <- dat_row[, c("state", "Date", "mean", "loci", "hici")]
    dat_row_all <- rbind.data.frame(dat_row_all, dat_row)
  }
  dat_row_all
}

median_CI_replicates_ratio <- function(dat) {
  
  dat_row_all <- NULL
  for (each_row in 1:nrow(dat)) {
    
    dat_row <- dat[each_row, ]
    
    if(any(!is.na(dat_row[, 3:(1+replicates)]))) {
      dat_row$median <- rowMedians(as.matrix(dat_row[, 3:(2+replicates)]), na.rm = TRUE)
      dat_row$loci <- quantile(dat_row[, 3:(2+replicates)], prob = 0.025, na.rm = TRUE)
      dat_row$hici <- quantile(dat_row[, 3:(2+replicates)], prob = 0.975, na.rm = TRUE) 
    } else {
      dat_row$median <- NA
      dat_row$loci <- NA
      dat_row$hici <- NA 
    }
    dat_row <- dat_row[, c("state", "Date", "median", "loci", "hici")]
    dat_row_all <- rbind.data.frame(dat_row_all, dat_row)
  }
  dat_row_all
}

rt_fun_replicate_parallel <- function(df, window_size, mcmc_control, seed) {
  
  # restrict to at least 12 infections to prevent to avert issues w/ posterior CV (see original EpiEstim paper)
  df <- as.matrix(df)
  date_names <- colnames(df)
  
  cumsum_inf <- rowCumsums(df)
  min12 <- min(which(cumsum_inf >= 12))
  df2 <- df[min12:length(df)]
  df2_dates <- date_names[min12:length(date_names)]
  
  # # 1st day with infections of variant to start the R estimate otherise R estimate artificially high
  # non0 <- min(which(df > 0)) 
  # 
  # # df filtered where there is the first case of variant to end of dataset
  # df2 <- df[non0:length(df)] 
  
  # if have negative infections, sets to 0 (can be caused by rolling average)
  df2[which(df2 < 0)] <- 0 
  
  # specify estimation window
  t_start <- seq(2, length(df2)-window_size) 
  t_end <- t_start + window_size
  
  # specify config statements
  config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                             std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                             n1=500,n2=50,t_start=t_start, t_end=t_end,
                             mcmc_control = mcmc_control, # for reproducibility
                             seed = seed))
  # n1 = number of pairs of mean and sd of the SI that are drawn 
  # n2 = size of the posterior sample drawn for each pair of mean, sd of SI
  
  # for uncertain si:
  # mean_si = positive real giving the average mean serial interval
  # std_mean_si = standard deviation of the distribution from which mean serial intervals are drawn
  # min_mean_si = lower bound of the distribution from which mean serial intervals are drawn
  # max_mean_si =  upper bound of the distribution from which mean serial intervals are drawn
  # std_si = non negative real giving the average standard deviation of the serial interval
  # std_std_si = standard deviation of the distribution from which standard deviations of the serial interval are drawn 
  # min_std_si = lower bound of the distribution from which standard deviations of the serial interval are drawn 
  # max_std_si = upper bound of the distribution from which standard deviations of the serial interval are drawn 
  
  # daily infection incidence vector
  df2 <- data.frame(I = df2)
  epiestim_output <- estimate_R(df2,
                                method="uncertain_si",
                                config = config)
  
  if (any(epiestim_output$si_distr < 0)) {
    print("warning: negative serial interval")
  }
  
  Rt_estimates <- data.frame(rt = epiestim_output$R$`Mean(R)`)
  
  # t_start vector must start after the 1st day of >12 inf incidence
  start_date <- as.Date(df2_dates[1]) + 1 + window_size
  Rt_estimates$Date <- start_date + 0:(length(Rt_estimates$rt) - 1)
  
  # join to full date range so can cbind dataframes
  df_date_names <- data.frame(Date = as.Date(date_names))
  
  Rt_estimates_all_dates <- df_date_names %>%
    dplyr::left_join(Rt_estimates, by = c("Date"))
  
  Rt_estimates_all_dates
}
