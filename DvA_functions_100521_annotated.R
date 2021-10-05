
# process initial GISAID dataset
readNextMeta_states <- function(next_meta_path){
  SS <- read_tsv(next_meta_path)
  SS <- SS[,c("Virus name", "Collection date", "Pango lineage", "AA Substitutions","Clade", "Location")]
  names(SS) <- c("seqName", "Date", "pango_lineage", "AAsubstitutions", "Clade", "Location")
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
check_outliers <- function(variant_df_vec, num_zero_days) {
  
  prob_log_list <- list()
  
  interest_col <- "I"
  
  for (state in states) {
    
    dat <- variant_df_vec[[state]]
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

rt_fun_12inf = function(df, name, window_size, mcmc_control, seed) {
  
  var_I <- "I"

  # restrict to at least 12 infections to prevent to avert issues w/ posterior CV (see original EpiEstim paper)
  df$Cumul_I <- rep(NA, nrow(df))
  df$Cumul_I <- cumsum(df[, var_I])
  min12 <- min(which(df[, "Cumul_I"] >= 12))
  date_min12 <- df[min12, "Date"]
  
  # df filtered where there are 12 cumulative variant infections to end of data
  df2 = df[min12:nrow(df), ] 
  
  # if have negative infections, sets to 0 (can be caused by rolling average)
  df2[which(df2[,var_I] < 0), var_I] <- 0 
  
  # specify estimation window
  t_start <- seq(2, nrow(df2)-window_size) 
  t_end <- t_start + window_size
  
  # specify config statements
  config <- make_config(list(mean_si = 5.2, std_mean_si = 1,min_mean_si = 2.2, max_mean_si = 8.2,
                             std_si = 4, std_std_si = 0.5,min_std_si = 2.5, max_std_si = 5.5,
                             n1=500,n2=50,t_start=t_start, t_end=t_end,
                             mcmc_control = mcmc_control, # for reproducibility
                             seed = seed)
                        # n1 = number of pairs of mean and sd of the SI that are drawn 
                        # n2 = size of the posterior sample drawn for each pair of mean, sd of SI
  )
  
  # daily infection incidence vector
  epiestim_output = estimate_R(df2[ , var_I],
                               method="uncertain_si",
                               config = config)
  
  Rt_estimates <- epiestim_output$R
  I_estimates <- data.frame(epiestim_output$I)
  I_estimates$Date <- df2$Date
  
  # adds back in days that were filtered out to match the days in the main dataframe
  start_date <- date_min12 + 1 + window_size # t_start vector must start after the 1st day of non-null incidence
  Rt_estimates$Date <- start_date + 0:(nrow(Rt_estimates) - 1)
  Rt_estimates <- left_join(Rt_estimates, I_estimates, by = c("Date"))
  
  Rt_estimates <- Rt_estimates %>%
    dplyr::rename(rt = `Mean(R)`,
                  rtlowci = `Quantile.0.25(R)`,
                  rtupci = `Quantile.0.975(R)`) %>%
    dplyr::select(Date, t_start, t_end, rt, rtlowci, rtupci)
  
  df_sub <- df[, c("state", "Date", "I")]
  
  merge <- df_sub %>% 
    dplyr::left_join(Rt_estimates, by = "Date") 
  
  merge
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
