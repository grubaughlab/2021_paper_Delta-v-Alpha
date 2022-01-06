#############################################################################################################################
# ABOUT ----------------------------------------------------------------------------------------------------
#############################################################################################################################

# This code generates plots related to sequencing coverage, variant frequencies, and logistic growth rates. 
## Figures 1-2 in the main text
## Figures S1 and S4 in the supplementary material

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
path <- getwd()
file_name <- "/DvA_functions.R"
source(full_path, echo=TRUE)

# Load GISAID metadata file 
metadata_date <- "2021-08-13"
gisaid_date <- paste("/", metadata_date, "/", sep="")
# next_meta <- readNextMeta_states(paste0(path, "metadata_", metadata_date, ".tsv.gz", sep="")) # run if 1st time running code
# save(next_meta, file = "next_meta.Rda") # save so don't have to re-process next time
path_data <- paste0(path, "data/")
load(paste0(path_data, "next_meta.Rda"))

# set the number of days for the emergence period
emergence_pd <- 90

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
# Percent Fully Vaccinated Over Time ----------------------------------------------
#############################################################################################################################

file_name <- "covidestim_state_2021-11-22.csv"
infect_import <- read.csv(paste0(path_data, file_name))

infect = infect_import %>%
  dplyr::select(state, date, infections) %>% # Median estimated infections/day, by date of infection
  # dplyr::select(state, date, infections, infections.lo, infections.hi) %>% # Median and 95% interval around estimated infections/day, by date of infection
  dplyr::filter(state %in% states, date >= "2021-01-01", date <= "2021-07-31") %>% # matching dates for var frequencies (removing dates w/ NAs at start & end caused by 7 day rolling average)
  dplyr::rename(Date = date) %>%
  dplyr::mutate(Date = as.Date(Date))

# get vx prop by state over time
# Source = https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc 
dat_vx <- read_csv(paste0(path_data, "COVID-19_Vaccinations_in_the_United_States_Jurisdiction.csv"))

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
vx_inf <- vx_inf + geom_text(x = as.Date("2021-02-15"), y = 235, color = "black", size = 2.5, aes(label = label), data = vx_inf_labels)
vx_inf_labels <- data.frame(state = c("Connecticut"), 
                            label = c("Delta Emergence Period"))
vx_inf <- vx_inf + geom_text(x = as.Date("2021-06-03"), y = 235, color = "black", size = 2.5, aes(label = label), data = vx_inf_labels)
vx_inf

#############################################################################################################################
# Main Figure 1: Sequencing Coverage & Variant Frequencies ----------------------------------------------
#############################################################################################################################
library(gridExtra)
library(ggpubr)

p <- ggarrange(perc_cases_seq_roll7 + guides(color = guide_legend(nrow=2,byrow=TRUE,override.aes = list(alpha = 0)),
                                             fill=guide_legend(nrow=2,byrow=TRUE,override.aes = list(alpha = 0))),
               freq_time + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
               labels = c("A", "B"),
               nrow = 1,
               widths = c(1.5, 1.5), 
               heights = c(5,5))

ggsave(paste(path_figures, "Fig_1.pdf", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Main Figure 2: Vx/Inf, Growth Rates, Slope  ----------------------------------------------
#############################################################################################################################

mylegend <- g_legend(growth_plots_emerge)
p <- grid.arrange(ggarrange(vx_inf + theme(legend.position="none"),
                            growth_plots_emerge + theme(legend.position="none"),
                            coef_bar_emerge + theme(legend.position="none"),
                            nrow = 1,
                            widths = c(1.3, 1, 0.85),
                            labels = c("A", "B", "C")),
                  mylegend,
                  nrow = 2,
                  heights = c(9.5, 0.5))

ggsave(paste(path_figures, "Fig_2.pdf", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Supplementary Figure S1: Vx/Inf/Date of Detection  ----------------------------------------------
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
  geom_text(size = 3, nudge_x = 2.5, nudge_y = -0.5) +
  xlab("Date of 1st Delta Detection") +
  xlim(as.Date("2021-03-15"), as.Date("2021-05-01")) +
  ylab("Percent of Population Fully Vaccinated at 1st Delta Detection")

p <- grid.arrange(ggarrange(emerge_gr_vx + theme(legend.position="none"),
                            emerge_gr_inf + theme(legend.position="none") ,
                            emerge_vx_date + theme(legend.position="none"),
                            labels = c("A", "B", "C")))

ggsave(paste(path_figures, "Fig_S1.pdf", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Supplementary Figure S4: Multiplicative Increase in R ----------------------------------------------
#############################################################################################################################

p <- fold_increase_emerge
ggsave(paste(path_figures, "Fig_S4.pdf", sep=""), p, height = 5, width =10)

#############################################################################################################################
# Summary Stats  ----------------------------------------------
#############################################################################################################################

# Sequencing Coverage ***********************************************************************

## maximum daily coverage by state
max_seq_cov_state <- ddply(cases_seq_date_7_all, .(state), summarize, Max_Seq = round(max(Percent_Seq7, na.rm = TRUE), digits = 2))
max_seq_cov_state[order(max_seq_cov_state$Max_Seq, decreasing = TRUE), ]

## date all states hit at least 5% seq coverage consistently
above_5p_seq_cov <- cases_seq_date_7_all
above_5p_seq_cov$Above_5 <- rep(NA, nrow(above_5p_seq_cov))
above_5p_seq_cov[which(above_5p_seq_cov$Percent_Seq7 >= 0.05), "Above_5"] <- "Y"
above_5p_seq_cov[which(above_5p_seq_cov$Percent_Seq7 < 0.05), "Above_5"] <- "N"

above_5p_seq_cov_first_time <- ddply(above_5p_seq_cov %>% dplyr::filter(Above_5 == "Y"), .(state), summarize, Min_5perc = min(Date, na.rm = TRUE))
above_5p_seq_cov <- above_5p_seq_cov %>% filter(Date >= max(above_5p_seq_cov_first_time$Min_5perc))
any(above_5p_seq_cov$Percent_Seq7 < 0.05) # check that no states dipped back below 0.05
check_5perc <- above_5p_seq_cov %>% dplyr::filter(Above_5 == "N", Date < "2021-07-15") # restrict due to delays in sequencing / reporting sequencing data

check_5perc$epidate <- epiweek(check_5perc$Date)
check_5perc <- check_5perc %>% dplyr::filter(Date > "2021-04-08")
ddply(check_5perc, .(state), summarize, Max_Date_Below_5p = max(Date))
ddply(check_5perc, .(state), nrow) # number of days below 5%
ddply(check_5perc, .(state), summarize, Mean_Cov = round(mean(Percent_Seq7), digits = 2)) # avg coverage in <5% days

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
df_delta <- df %>% dplyr::filter(variant_category == "Delta", epidate == "2021-07-26") # last week for which have data for all states due to  seq / reporting delay
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

# Multiplicative Rt Increase  ***********************************************************************
mean_fold_inc <- ddply(model_df, .(lineage), summarize, Mean_Trans = round(mean(transmissibility), digits = 2))
mean_fold_inc
mean_fold_inc_state <- ddply(model_df, .(state, lineage), summarize, Mean_Trans = round(mean(transmissibility), digits = 2))
mean_fold_inc_state_cast <- reshape2::dcast(mean_fold_inc_state, state ~ lineage, value.var = "Mean_Trans")
mean_fold_inc_state_cast[order(mean_fold_inc_state_cast$Delta, decreasing = TRUE), ]

