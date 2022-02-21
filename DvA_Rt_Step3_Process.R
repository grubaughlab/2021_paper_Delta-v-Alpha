#############################################################################################################################
# ABOUT ----------------------------------------------------------------------------------------------------
#############################################################################################################################

# This code takes the truncated raw and summarized Rt estimate data and prepares it for plotting.
## Figure 3 in the main text
## Figures S5-6 in the supplementary text

#############################################################################################################################
# LOAD LIBRARIES & FILES ----------------------------------------------------------------------------------------------------
#############################################################################################################################

library(RColorBrewer)
library(ggplot2)
library(scales)

# number of replicates in files
replicates <- 1000

# load functions
path <- getwd()
file_name <- "/DvA_functions.R"
full_path <- paste0(path, file_name)
source(full_path, echo=TRUE)

# path for where to save figures
path_figures <- paste0(path, "figures/")

# path from where to grab data files
path_data <- paste0(path, "data/")

# specify selected states
states <- c("Connecticut", "Massachusetts", "Maine", "New Hampshire", "Rhode Island", "Vermont")

# select colors for plotting
all_dark2_colors <- brewer.pal(8, "Dark2")
customPalette <- c(all_dark2_colors[1], all_dark2_colors[3], all_dark2_colors[7]) # Alpha = green; Delta = purple; Other = brown

#############################################################################################################################
# Mean Rt Over Time ----------------------------------------------
#############################################################################################################################

# load Rt estimate replicates
load(paste0(path_data, "all_rt_trunc_all_", replicates, ".rda"))

# plot Rt estimates mean and CI - all variant categories
max_plot <- ceiling(max(all_rt_trunc_all$mean, na.rm = TRUE))
rt_plot_all <- ggplot(all_rt_trunc_all, aes(Date, mean)) +
  geom_line(aes(color = Variant_Category)) +
  geom_ribbon(aes(ymin=loci, ymax=hici, fill = Variant_Category), linetype = 0, alpha = 0.2) +
  scale_fill_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) +
  scale_color_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) + 
  scale_y_continuous(breaks = seq(0, max_plot, by = 1), limits = c(0, max_plot)) +
  facet_wrap(~state, ncol = 1) +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Rt") +
  theme(legend.position = "bottom", legend.box = "horizontal")  +
  labs(color = "Variant Category", fill = "Variant Category") 
rt_plot_all

# full CI
rt_plot_all_CI <- ggplot(all_rt_trunc_all, aes(Date, mean)) +
  geom_line(aes(color = Variant_Category)) +
  geom_ribbon(aes(ymin=loci, ymax=hici, fill = Variant_Category), linetype = 0, alpha = 0.2) +
  scale_fill_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) +
  scale_color_manual(values = customPalette, labels = c("Alpha", "Delta", "Other")) + 
  facet_wrap(~state, ncol = 1, scales = "free_y") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Rt") +
  theme(legend.position = "bottom", legend.box = "horizontal")  +
  labs(color = "Variant Category", fill = "Variant Category") 
rt_plot_all_CI

#############################################################################################################################
# Mean Rt Ratio Over Time ----------------------------------------------
#############################################################################################################################

# load Rt estimate replicates
load(paste0(path_data, "all_rt_trunc_all_raw_", replicates, ".rda"))

alpha_rt <- all_rt_trunc_all_raw %>% dplyr::filter(Variant_Category == "Alpha")
delta_rt <- all_rt_trunc_all_raw %>% dplyr::filter(Variant_Category == "Delta")
other_rt <- all_rt_trunc_all_raw %>% dplyr::filter(Variant_Category == "Other")

# calculate Rt ratio - Delta:Alpha 
all_ratios <- NULL
for (each_rep in 1:replicates) {
  
  alpha_rt_rep <- alpha_rt[, c("Date", "state", paste0("rt_", each_rep))]
  delta_rt_rep <- delta_rt[, c("Date", "state", paste0("rt_", each_rep))]
  colnames(alpha_rt_rep) <- c("Date", "state", "rt_alpha")
  colnames(delta_rt_rep) <- c("Date", "state", "rt_delta")
  
  delta_to_alpha <- left_join(delta_rt_rep, alpha_rt_rep, by = c("state", "Date"))
  delta_to_alpha$DvA_ratio <-  delta_to_alpha$rt_delta / delta_to_alpha$rt_alpha
  delta_to_alpha <- delta_to_alpha[, c("Date", "state", "DvA_ratio")]
  colnames(delta_to_alpha) <- c("Date", "state", paste0("DvA_ratio_", each_rep))
  if (each_rep == 1) {
    all_ratios <- delta_to_alpha
  } else {
    all_ratios <- left_join(all_ratios, delta_to_alpha, by = c("Date" = "Date", "state" = "state"))
  }
  all_ratios
}

# calculate mean and 95% CI 
library(matrixStats)
# all_ratios_summary <- median_CI_replicates_ratio(all_ratios) # if opted to use median instead of mean
all_ratios_summary <- mean_CI_replicates_ratio(all_ratios)
all_ratios_summary$Comparison <- "Delta:Alpha"
max_level <- max(all_ratios_summary$mean, na.rm = TRUE)

rt_ratio_ad <- all_ratios_summary %>%
  ggplot(aes(x = Date, y = mean, color = Comparison, fill = Comparison)) +
  geom_line() +
  facet_wrap(~state, ncol = 1) +
  scale_fill_manual(values = "dark grey") +
  scale_color_manual(values = "dark grey") +
  geom_ribbon(aes(ymin=loci, ymax=hici), linetype = 0, alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Mean Rt Ratio") +
  scale_y_continuous(breaks = seq(0, max_level, by = 2), limits = c(0, max_level)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")  +
  theme(legend.position = "bottom", legend.box = "horizontal") 

rt_ratio_ad

# full CI
rt_ratio_ad_CI <- all_ratios_summary %>%
  ggplot(aes(x = Date, y = mean, color = Comparison, fill = Comparison)) +
  geom_line() +
  facet_wrap(~state, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = "dark grey") +
  scale_color_manual(values = "dark grey") +
  geom_ribbon(aes(ymin=loci, ymax=hici), linetype = 0, alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Mean Rt Ratio") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_x_date(date_breaks = "1 month", date_labels = "%b")  +
  theme(legend.position = "bottom", legend.box = "horizontal") 

rt_ratio_ad_CI

#############################################################################################################################
# Rt Replicates ----------------------------------------------
#############################################################################################################################

## format for plotting
selected_columns <- rep(paste0("rt_", seq(1:replicates)))
alpha_rt_melt <- reshape2::melt(alpha_rt, id = c("Date", "state"), measure.vars = selected_columns)
delta_rt_melt <- reshape2::melt(delta_rt, id = c("Date", "state"), measure.vars = selected_columns)
other_rt_melt <- reshape2::melt(other_rt, id = c("Date", "state"), measure.vars = selected_columns)

rt_plot_alpha <- ggplot(alpha_rt_melt, aes(Date, value)) +
  geom_line(aes(color = variable), size = 0.2) +
  facet_wrap(~state, ncol = 1) +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Mean Rt (1000 samples)") +
  theme(legend.position = "none")

rt_plot_delta <- ggplot(delta_rt_melt, aes(Date, value)) +
  geom_line(aes(color = variable), size = 0.2) +
  facet_wrap(~state, ncol = 1, scales = "free_y") +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Mean Rt (1000 samples)") +
  theme(legend.position = "none")

rt_plot_other <- ggplot(other_rt_melt, aes(Date, value)) +
  geom_line(aes(color = variable), size = 0.2) +
  facet_wrap(~state, ncol = 1) +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = 1, color = "black", size = 0.3) +
  ylab("Estimated Mean Rt (1000 samples)") +
  theme(legend.position = "none")

# Add co-circulation period w/ estimated infections
load(paste0(path_data, "cocirc_plot_data.rda"))
p_cocirc <- ggplot(trunc_inf_all_pop_mean) + 
  geom_bar(aes(x = Date, y = Mean_Inf_per_100K, fill = Variant_Category, color = Variant_Category), stat = "identity") +
  scale_fill_manual(values = customPalette, name = "Variant Category") +
  scale_color_manual(values = customPalette, name = "Variant Category") +
  theme_bw() +
  theme(
    axis.title.y = element_text(color = "black"),
    axis.title.y.right = element_text(color = "black")) + 
  facet_wrap(~state, scales = "fixed", ncol = 1) +
  theme(legend.position="bottom") +
  ylab("Mean Estimated Infections per 100K Population") 

#############################################################################################################################
# Plots ----------------------------------------------
#############################################################################################################################

library(ggpubr)

# plot mean Rt and ratios
p <- ggarrange(rt_plot_all +guides(fill=guide_legend(nrow=2,byrow=TRUE)), 
               rt_ratio_ad + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
               labels = c("A", "B"),
               ncol = 2,
               nrow = 1,
               widths = c(1.5, 1.5), 
               heights = c(5, 5))

ggsave(paste(path_figures, "Fig_3.pdf", sep=""), p, height = 10, width =10)

# plot mean Rt and ratios - full CI
p <- ggarrange(p_cocirc + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
               rt_plot_all_CI +guides(fill=guide_legend(nrow=2,byrow=TRUE)), 
               rt_ratio_ad_CI + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
               labels = c("A", "B", "C"),
               ncol = 3,
               nrow = 1,
               widths = c(1.5, 1.5, 1.5), 
               heights = c(5, 5, 5))

ggsave(paste(path_figures, "Fig_S3.pdf", sep=""), p, height = 10, width =10)

# plot raw runs by sample
# p <- ggarrange(rt_plot_alpha +guides(fill=guide_legend(nrow=2,byrow=TRUE)), 
#                rt_plot_delta + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
#                rt_plot_other + guides(fill=guide_legend(nrow=2,byrow=TRUE)),
#                labels = c("A", "B", "C"),
#                ncol = 3,
#                nrow = 1,
#                widths = c(1.5, 1.5, 1.5), 
#                heights = c(5, 5, 5))
# 
# ggsave(paste(path_figures, "Fig_S6.pdf", sep=""), p, height = 10, width =10)

#############################################################################################################################
# Summary Statistics  ----------------------------------------------
#############################################################################################################################

## Rt 
mean_rt_full_period <- ddply(all_rt_trunc_all, .(state, Variant_Category), summarize, Mean_Rt = round(mean(mean, na.rm = TRUE), digits = 2))
mean_rt_full_period_cast <- reshape2::dcast(mean_rt_full_period, state ~ Variant_Category, value.var = "Mean_Rt")
mean_rt_full_period_cast[order(mean_rt_full_period_cast$Other, decreasing = TRUE), ]
# Other - 0.87 (Vermont, Rhode Island, Massachusetts) to 0.91 (Maine)

mean_rt_full_period_delta <- ddply(all_rt_trunc_all, .(Variant_Category), summarize, Mean_Rt = round(mean(mean, na.rm = TRUE), digits = 2))
mean_rt_full_period_delta %>% dplyr::filter(Variant_Category == "Delta") # 1.4
mean_rt_full_period_cast[order(mean_rt_full_period_cast$Delta, decreasing = TRUE), ] # 1.27 (New Hampshire) to 1.65 (Vermont)

## Alpha mean Rt pre- and post-Delta
load(paste0(path_data, "cocirc_dates.rda"))

dat_alpha <- all_rt_trunc_all %>% dplyr::filter(Variant_Category == "Alpha")
dat_alpha$Period <- NA

dat_alpha_all <- list()
for (each_state in unique(dat_alpha$state)) {
  
  dat_state <- dat_alpha %>% dplyr::filter(state == each_state)
  combo_dates_state <- combo_dates  %>% dplyr::filter(state == each_state)
  dat_state[which(dat_state$Date < combo_dates_state$Delta_Start), "Period"] <-  "Pre_Delta"
  dat_state[which(dat_state$Date >= combo_dates_state$Delta_Start), "Period"] <-  "Post_Delta"
  dat_alpha_all[[each_state]] <- dat_state
}

dat_alpha_all <- do.call(rbind.data.frame, dat_alpha_all)

mean_rt_pre_post_delta <- ddply(dat_alpha_all, .(Variant_Category, Period), summarize, Mean_Rt = round(mean(mean, na.rm = TRUE), digits = 2))
mean_rt_pre_post_delta

# Rt ratio summary

## mean overall for each state
mean_Rt_ratio <- ddply(all_ratios_summary, .(state, Comparison), summarize, Mean_Rt_Ratio = round(mean(mean, na.rm = TRUE), digits = 2))
mean_Rt_ratio_cast <- reshape2::dcast(mean_Rt_ratio, state ~ Comparison, value.var = "Mean_Rt_Ratio")
mean_Rt_ratio_cast[order(mean_Rt_ratio_cast$`Delta:Alpha`, decreasing = TRUE), ]

# 1.63 - 2.67

## range for each state
min_max_Rt_ratio <- ddply(all_ratios_summary, .(state, Comparison), summarize, 
                          Min_Ratio = round(min(mean, na.rm = TRUE), digits = 2), 
                          Max_Ratio = round(max(mean, na.rm = TRUE), digits = 2))

min_max_Rt_ratio[order(min_max_Rt_ratio$Min_Ratio, decreasing = TRUE), ]








