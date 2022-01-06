#############################################################################################################################
# ABOUT ----------------------------------------------------------------------------------------------------
#############################################################################################################################

# This code analyzes and plots the PCR CT value data.
## Figure 4 in the main text
## Figures S7-10 in the supplementary text

#############################################################################################################################
# LOAD LIBRARIES & FILES ----------------------------------------------------------------------------------------------------
#############################################################################################################################

library(tidyverse)
library(dplyr)
library(ggpubr)
library(RColorBrewer)
library(readxl)

# set colors for plots
numGroups <- 4
customPalette <- brewer.pal(numGroups, "Dark2")
customPalette <- c(customPalette[1], customPalette[3], customPalette[4], customPalette[2])

path <- getwd()
path_data <- paste0(path, "data/CT/")
path_figures <- paste0(path, "figures/")

# read in data
CT_Yale <- read_csv(paste0(path_data, "Yale_Clinical_Virology/GLab_SC2_sequencing_data_082421.csv"))
CT_Jax <- read_xlsx(paste0(path_data, "CT_JAX/Delta vs Alpha - Ct data to send.xlsx"))
MA_MGB <- read_xlsx(paste0(path_data, "Mass_General_Brigham/Final_ct_values_alpha_delta.xlsx"))
ME_HETL_Jax <- read_csv(paste0(path_data, "ME_HETL_JAX/Ct_information_02092021.csv"))

# set variant categories
alpha_lineages <- c("B.1.1.7")
delta_lineages <- c("B.1.617.2", paste0("AY.", seq(from = 1, to = 25, by = 1)), "AY.3.1")

# set analysis time windows
alpha_start <- "2021-01-01"
alpha_end <- "2021-03-31"
delta_start <- "2021-06-01"
delta_end <- "2021-08-31"

# Connecticut / Yale --------------------------------------------------------------

# format data
CT_Yale$`Collection-date` <- as.Date(CT_Yale$`Collection-date`)
CT_Yale$`Yale-N1(FAM)` <- as.numeric(CT_Yale$`Yale-N1(FAM)`)

CT_Yale <- CT_Yale %>% dplyr::rename(Collection_Date = `Collection-date`,
                                     CT = `Yale-N1(FAM)`,
                                     State = `Division (state)`)

# subset data
CT_Yale_sub <- CT_Yale %>% 
  dplyr::filter(Source == "Yale Clinical Virology Lab", 
         Collection_Date >= alpha_start, 
         Collection_Date <= delta_end,
         State == "Connecticut",
         Lineage %in% c(alpha_lineages, delta_lineages))  %>% 
  drop_na(CT) %>% 
  dplyr::select(Source, Collection_Date, State, CT, Lineage)

# specify Alpha v Delta
dat <- CT_Yale_sub
dat$Variant_Category <- rep(NA, nrow(dat))
dat[which(dat$Lineage %in% alpha_lineages), "Variant_Category"] <- "Alpha"
dat[which(dat$Lineage %in% delta_lineages), "Variant_Category"] <- "Delta"

alpha_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Alpha", Collection_Date >= alpha_start, Collection_Date <= alpha_end) 

alpha_CT <- alpha_CT %>% 
  dplyr::filter(Collection_Date >= "2021-01-27") # when Yale started using current assay; drops 4 sequences

delta_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Delta", Collection_Date >= delta_start, Collection_Date <= delta_end) 

all_CT <- rbind.data.frame(alpha_CT, delta_CT)

# generate summary stats
ss_alpha <- paste0("N=", nrow(all_CT %>% filter(Variant_Category == "Alpha")))
ss_delta <- paste0("N=", nrow(all_CT %>% filter(Variant_Category == "Delta")))

alpha_CT <- all_CT %>% dplyr::filter(Variant_Category == "Alpha")
delta_CT <- all_CT %>% dplyr::filter(Variant_Category == "Delta")

mean_alpha <- paste0("Mean=", round(mean(alpha_CT$CT), digits = 2))
mean_delta <- paste0("Mean=", round(mean(delta_CT$CT), digits = 2))

median_alpha <- paste0("Median=", round(median(alpha_CT$CT), digits = 2))
median_delta <- paste0("Median=", round(median(delta_CT$CT), digits = 2))

summary_alpha <- paste0("Alpha (Jan-Mar) ", ss_alpha, " ", mean_alpha, " ", median_alpha)
summary_delta <- paste0("Delta (Jun-Aug) ", ss_delta, " ", mean_delta, " ", median_delta)

# plot
addline_format <- function(x,...){
  gsub('\\s','\n',x)
}

min(all_CT$CT)
max(all_CT$CT)

p  <- ggplot(all_CT, aes(x=Variant_Category, y=CT, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("PCR CT Value") +
  ggtitle("Yale University (Connecticut)") +
  theme_bw() +
  scale_x_discrete(labels = addline_format(c(summary_alpha, summary_delta))) + 
  scale_y_continuous(trans = "reverse", breaks = seq(10, 40, by = 10), limits = c(43, 8)) +
  theme(axis.title.x=element_blank(), legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) +
  stat_compare_means(data = all_CT, mapping = aes(x= Variant_Category, y = CT), 
                     comparisons = list(c("Alpha", "Delta")), method = "t.test", label = "p.signif", vjust = 0.25, hide.ns = FALSE) 

p
plot_CT_Yale <- p

p  <- ggplot(all_CT, aes(x=CT, fill = Variant_Category)) +
  geom_histogram(aes(fill = Variant_Category, color = Variant_Category), alpha = 0.9, position = "identity") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("Number of Samples") +
  xlab("PCR CT Value") +
  ggtitle("Yale University (Connecticut)") +
  theme_bw() +
  labs(color = "Variant Category", fill = "Variant Category")
p
plot_CT_Yale_hist <- p

mu <- ddply(all_CT, .(Variant_Category), summarize, grp.mean=mean(CT))
p  <- ggplot(all_CT, aes(x=CT, fill = Variant_Category)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("Density") +
  xlab("PCR CT Value") +
  ggtitle("Yale University (Connecticut)") +
  theme_bw() +
  labs(color = "Variant Category", fill = "Variant Category") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Variant_Category),
             linetype="dashed")
p
plot_CT_Yale_dens <- p

# #### Supplementary Figure S3: Alpha & Delta Values over Time -------------------------------

all_CT_alpha <- dat %>% dplyr::filter(Variant_Category == "Alpha")
all_CT_alpha$Month <- lubridate::month(all_CT_alpha$Collection_Date, label = TRUE)

p <- ggplot(all_CT_alpha, aes(x=Month, y=CT, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("PCR CT Value") +
  xlab("Month") +
  ggtitle("Yale University (Connecticut) - Alpha") +
  theme_bw() +
  scale_y_continuous(trans = "reverse") +
  theme(legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) 

p
plot_CT_Yale_Alpha_Time <- p

anova_alpha <- aov(CT ~ as.factor(Month), data = all_CT_alpha)
summary(anova_alpha)
# plot(anova_alpha, 1)
# plot(anova_alpha, 2)

tukey_alpha <- TukeyHSD(anova_alpha)
df_hsd_alpha <- data.frame(tukey_alpha$`as.factor(Month)`)
df_hsd_alpha[order(df_hsd_alpha$p.adj, decreasing = FALSE), ]

all_CT_delta <- dat %>% dplyr::filter(Variant_Category == "Delta")
all_CT_delta$Month <- lubridate::month(all_CT_delta$Collection_Date, label = TRUE)

p <- ggplot(all_CT_delta, aes(x=Month, y=CT, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), delta = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("PCR CT Value") +
  xlab("Month") +
  ggtitle("Yale University (Connecticut) - Delta") +
  theme_bw() +
  scale_y_continuous(trans = "reverse") +
  theme(legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) 

p
plot_CT_Yale_Delta_Time <- p

anova_delta <- aov(CT ~ as.factor(Month), data = all_CT_delta)
summary(anova_delta)
# plot(anova_delta, 1)
# plot(anova_delta, 2)

# tukey_delta <- TukeyHSD(anova_delta)
# df_hsd_delta <- data.frame(tukey_delta$`as.factor(Month)`)
# df_hsd_delta[order(df_hsd_delta$p.adj, decreasing = FALSE), ]

# #### Supplementary Figure S4: Viral Loads ------------------------------------------

slope <- -3.762190476
intercept <- 43.3754127

all_CT <- all_CT %>%
  dplyr::mutate(log10_ge_ul = (CT - intercept) / slope, 
                ge_ul = 10^log10_ge_ul,
                elution_ul = rep(75, n()),
                extraction_ul = rep(300, n()),
                ge_ml = (ge_ul * elution_ul) * (1000 / extraction_ul))

all_CT$log10_ge_ml <- log10(all_CT$ge_ml) # NOTE: ml, not ul

# generate summary stats
ss_alpha <- paste0("N=", nrow(all_CT %>% filter(Variant_Category == "Alpha")))
ss_delta <- paste0("N=", nrow(all_CT %>% filter(Variant_Category == "Delta")))

alpha_CT <- all_CT %>% dplyr::filter(Variant_Category == "Alpha")
delta_CT <- all_CT %>% dplyr::filter(Variant_Category == "Delta")

mean_alpha <- paste0("Mean=", round(mean(alpha_CT$log10_ge_ml), digits = 2))
mean_delta <- paste0("Mean=", round(mean(delta_CT$log10_ge_ml), digits = 2))

mean_alpha_nonlog <- round(mean(alpha_CT$ge_ml), digits = 2)
mean_delta_nonlog <- round(mean(delta_CT$ge_ml), digits = 2)
round((mean_delta_nonlog/mean_alpha_nonlog), digits = 2)
  
median_alpha <- paste0("Median=", round(median(alpha_CT$log10_ge_ml), digits = 2))
median_delta <- paste0("Median=", round(median(delta_CT$log10_ge_ml), digits = 2))

summary_alpha <- paste0("Alpha (Jan-Mar) ", ss_alpha, " ", mean_alpha, " ", median_alpha)
summary_delta <- paste0("Delta (Jun-Aug) ", ss_delta, " ", mean_delta, " ", median_delta)

p  <- ggplot(all_CT, aes(x=Variant_Category, y=log10_ge_ml, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("Log10 Virus RNA Copies per mL") +
  ggtitle("Yale University (Connecticut)") +
  theme_bw() +
  scale_x_discrete(labels = addline_format(c(summary_alpha, summary_delta))) + 
  theme(axis.title.x=element_blank(), legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) +
  stat_compare_means(data = all_CT, mapping = aes(x= Variant_Category, y = log10_ge_ml), 
                     comparisons = list(c("Alpha", "Delta")), method = "t.test", label = "p.signif", vjust = 0.25, hide.ns = FALSE) 

p
plot_CT_Yale_Std <- p

p <- grid.arrange(ggarrange(plot_CT_Yale_Std))
ggsave(paste(path_figures, "Fig_S8.pdf", sep=""), p, height = 5, width =5)

# Connecticut / JAX Data ----------------------------------------------------------------

# format data
CT_Jax$`Collection Date` <- as.Date(CT_Jax$`Collection Date`)
CT_Jax <- CT_Jax %>% dplyr::rename(Collection_Date = `Collection Date`,
                                   Lineage = `Lineage Call`,
                                   CT = `Ct Value`)

# subset data - drop SGTF samples not sequenced
CT_Jax_sub <- CT_Jax %>% dplyr::filter(Class %in% c("Sequenced - B.1.1.7", "Sequenced - B.1.617.2"))

# specify Alpha v Delta
dat <- CT_Jax_sub
dat$Variant_Category <- rep(NA, nrow(dat))
dat[which(dat$Lineage %in% alpha_lineages), "Variant_Category"] <- "Alpha"
dat[which(dat$Lineage %in% delta_lineages), "Variant_Category"] <- "Delta"

alpha_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Alpha", Collection_Date >= alpha_start, Collection_Date <= alpha_end) 

delta_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Delta", Collection_Date >= delta_start, Collection_Date <= delta_end) 

all_CT <- rbind.data.frame(alpha_CT, delta_CT)

# generate summary stats
ss_alpha <- paste0("N=", nrow(all_CT %>% dplyr::filter(Variant_Category == "Alpha")))
ss_delta <- paste0("N=", nrow(all_CT %>% dplyr::filter(Variant_Category == "Delta")))

alpha_CT <- all_CT %>% dplyr::filter(Variant_Category == "Alpha")
delta_CT <- all_CT %>% dplyr::filter(Variant_Category == "Delta")

mean_alpha <- paste0("Mean=", round(mean(alpha_CT$CT), digits = 2))
mean_delta <- paste0("Mean=", round(mean(delta_CT$CT), digits = 2))

median_alpha <- paste0("Median=", round(median(alpha_CT$CT), digits = 2))
median_delta <- paste0("Median=", round(median(delta_CT$CT), digits = 2))

summary_alpha <- paste0("Alpha (Jan-Mar) ", ss_alpha, " ", mean_alpha, " ", median_alpha)
summary_delta <- paste0("Delta (Jun-Aug) ", ss_delta, " ", mean_delta, " ", median_delta)

min(all_CT$CT)
max(all_CT$CT)

p <- ggplot(all_CT, aes(x=Variant_Category, y=CT, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("PCR CT Value") +
  ggtitle("Jackson Laboratory (Connecticut)") +
  theme_bw() +
  scale_x_discrete(labels = addline_format(c(summary_alpha, summary_delta))) + 
  scale_y_continuous(trans = "reverse", breaks = seq(10, 40, by = 10), limits = c(43, 8)) +
  theme(axis.title.x=element_blank(), legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) +
  stat_compare_means(data = all_CT, mapping = aes(x= Variant_Category, y = CT), 
                     comparisons = list(c("Alpha", "Delta")), method = "t.test", label = "p.signif", vjust = 0.25, hide.ns = FALSE) 

p 
plot_CT_Jax <- p

mu <- ddply(all_CT, .(Variant_Category), summarize, grp.mean=mean(CT))
p  <- ggplot(all_CT, aes(x=CT, fill = Variant_Category)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("Density") +
  xlab("PCR CT Value") +
  ggtitle("Jackson Laboratory (Connecticut)") +
  theme_bw() +
  labs(color = "Variant Category", fill = "Variant Category") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Variant_Category),
             linetype="dashed")
p
plot_CT_Jax_dens <- p

# Massachusetts / MGB Data ----------------------------------------------------------------

# format data 
MA_MGB <- MA_MGB %>% dplyr::rename(Collection_Date = coll.date,
                                   Lineage = pango_lineage,
                                   CT_orf1 = orf1.ct,
                                   CT_e = e.ct, 
                                   Sample_Type = swab_source)
MA_MGB$Collection_Date <- as.Date(MA_MGB$Collection_Date, tryFormats = "%m/%d/%y")

MA_MGB <- MA_MGB %>% dplyr::filter(Sample_Type %in% c("NP", "NPSW", "NASL"))

# specify Alpha v Delta
dat <- MA_MGB
dat$Variant_Category <- rep(NA, nrow(dat))
dat[which(dat$Lineage %in% alpha_lineages), "Variant_Category"] <- "Alpha"
dat[which(dat$Lineage %in% delta_lineages), "Variant_Category"] <- "Delta"

alpha_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Alpha", Collection_Date >= alpha_start, Collection_Date <= alpha_end) 

delta_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Delta", Collection_Date >= delta_start, Collection_Date <= delta_end) 

all_CT <- rbind.data.frame(alpha_CT, delta_CT)

# generate summary stats
count(is.na(all_CT$CT_orf1)) # NA for several samples
count(is.na(all_CT$CT_e))

ss_alpha <- paste0("N=", nrow(all_CT %>% dplyr::filter(Variant_Category == "Alpha")))
ss_delta <- paste0("N=", nrow(all_CT %>% dplyr::filter(Variant_Category == "Delta")))

alpha_CT <- all_CT %>% dplyr::filter(Variant_Category == "Alpha")
delta_CT <- all_CT %>% dplyr::filter(Variant_Category == "Delta")

mean_alpha <- paste0("Mean=", round(mean(alpha_CT$CT_e), digits = 2))
mean_delta <- paste0("Mean=", round(mean(delta_CT$CT_e), digits = 2))

median_alpha <- paste0("Median=", round(median(alpha_CT$CT_e), digits = 2))
median_delta <- paste0("Median=", round(median(delta_CT$CT_e), digits = 2))

summary_alpha <- paste0("Alpha (Jan-Mar) ", ss_alpha, " ", mean_alpha, " ", median_alpha)
summary_delta <- paste0("Delta (Jun-Aug) ", ss_delta, " ", mean_delta, " ", median_delta)

min(all_CT$CT_e)
max(all_CT$CT_e)

p  <- ggplot(all_CT, aes(x=Variant_Category, y=CT_e, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("PCR CT Value") +
  ggtitle("Mass General Brigham (Massachusetts)") +
  theme_bw() +
  scale_x_discrete(labels = addline_format(c(summary_alpha, summary_delta))) + 
  scale_y_continuous(trans = "reverse", breaks = seq(10, 40, by = 10), limits = c(43, 8)) +
  theme(axis.title.x=element_blank(), legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) +
  stat_compare_means(data = all_CT, mapping = aes(x= Variant_Category, y = CT_e), 
                     comparisons = list(c("Alpha", "Delta")), method = "t.test", label = "p.signif", vjust = 0.25, hide.ns = FALSE) 

p 
plot_MA_MGB <- p

mu <- ddply(all_CT, .(Variant_Category), summarize, grp.mean=mean(CT_e))
p  <- ggplot(all_CT, aes(x=CT_e, fill = Variant_Category)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("Density") +
  xlab("PCR CT Value") +
  ggtitle("Mass General Brigham (Massachusetts)") +
  theme_bw() +
  labs(color = "Variant Category", fill = "Variant Category") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Variant_Category),
             linetype="dashed")
p
plot_MA_MGB_dens <- p

# #### Supplementary FigureS5: ORF1 Target ------------------------------------------
all_CT_orf1 <- all_CT %>% dplyr::filter(CT_orf1 != "NA") # NAs in orf1 target data
all_CT_orf1$CT_orf1 <- as.numeric(all_CT_orf1$CT_orf1)

ss_alpha <- paste0("N=", nrow(all_CT_orf1 %>% dplyr::filter(Variant_Category == "Alpha")))
ss_delta <- paste0("N=", nrow(all_CT_orf1 %>% dplyr::filter(Variant_Category == "Delta")))

alpha_CT <- all_CT_orf1 %>% dplyr::filter(Variant_Category == "Alpha")
delta_CT <- all_CT_orf1 %>% dplyr::filter(Variant_Category == "Delta")

mean_alpha <- paste0("Mean=", round(mean(alpha_CT$CT_orf1), digits = 2))
mean_delta <- paste0("Mean=", round(mean(delta_CT$CT_orf1), digits = 2))

median_alpha <- paste0("Median=", round(median(alpha_CT$CT_orf1), digits = 2))
median_delta <- paste0("Median=", round(median(delta_CT$CT_orf1), digits = 2))

summary_alpha <- paste0("Alpha (Jan-Mar) ", ss_alpha, " ", mean_alpha, " ", median_alpha)
summary_delta <- paste0("Delta (Jun-Aug) ", ss_delta, " ", mean_delta, " ", median_delta)

min(all_CT_orf1$CT_orf1)
max(all_CT_orf1$CT_orf1)

p  <- ggplot(all_CT_orf1, aes(x=Variant_Category, y=CT_orf1, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("PCR CT Value") +
  ggtitle("Mass General Brigham (Massachusetts)") +
  theme_bw() +
  scale_x_discrete(labels = addline_format(c(summary_alpha, summary_delta))) + 
  scale_y_continuous(trans = "reverse", breaks = seq(10, 40, by = 10), limits = c(43, 8)) +
  theme(axis.title.x=element_blank(), legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) +
  stat_compare_means(data = all_CT_orf1, mapping = aes(x= Variant_Category, y = CT_orf1), 
                     comparisons = list(c("Alpha", "Delta")), method = "t.test", label = "p.signif", vjust = 0.25, hide.ns = FALSE) 

p 
plot_MA_MGB_orf1 <- p
ggsave(paste(path_figures, "Fig_S9.pdf", sep=""), height = 5, width =5)

# Maine / HETL / Jax Data ----------------------------------------------------------------

# format data
ME_HETL_Jax$`Date of Collection` <- as.Date(ME_HETL_Jax$`Date of Collection`, tryFormats = "%m/%d/%y")
ME_HETL_Jax <- ME_HETL_Jax %>% dplyr::rename(Collection_Date = `Date of Collection`,
                                             Lineage = Pangolin,
                                             CT = `Ct (raw data)`)

# subset data - only using lab 1
ME_HETL_Jax_sub <- ME_HETL_Jax %>% dplyr::filter(Lab == 1)
ME_HETL_Jax_sub$CT <- as.numeric(ME_HETL_Jax_sub$CT) 

# specify Alpha v Delta
dat <- ME_HETL_Jax_sub
dat$Variant_Category <- rep(NA, nrow(dat))
dat[which(dat$Lineage %in% alpha_lineages), "Variant_Category"] <- "Alpha"
dat[which(dat$Lineage %in% delta_lineages), "Variant_Category"] <- "Delta"

alpha_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Alpha", Collection_Date >= alpha_start, Collection_Date <= alpha_end) 

delta_CT <- dat %>% 
  dplyr::filter(Variant_Category == "Delta", Collection_Date >= delta_start, Collection_Date <= delta_end) 

all_CT <- rbind.data.frame(alpha_CT, delta_CT)

# generate summary stats
ss_alpha <- paste0("N=", nrow(all_CT %>% dplyr::filter(Variant_Category == "Alpha")))
ss_delta <- paste0("N=", nrow(all_CT %>% dplyr::filter(Variant_Category == "Delta")))

alpha_CT <- all_CT %>% dplyr::filter(Variant_Category == "Alpha")
delta_CT <- all_CT %>% dplyr::filter(Variant_Category == "Delta")

mean_alpha <- paste0("Mean=", round(mean(alpha_CT$CT), digits = 2))
mean_delta <- paste0("Mean=", round(mean(delta_CT$CT), digits = 2))

median_alpha <- paste0("Median=", round(median(alpha_CT$CT), digits = 2))
median_delta <- paste0("Median=", round(median(delta_CT$CT), digits = 2))

summary_alpha <- paste0("Alpha (Jan-Mar) ", ss_alpha, " ", mean_alpha, " ", median_alpha)
summary_delta <- paste0("Delta (Jun-Aug) ", ss_delta, " ", mean_delta, " ", median_delta)

min(all_CT$CT)
max(all_CT$CT)

p  <- ggplot(all_CT, aes(x=Variant_Category, y=CT, color = Variant_Category)) +
  geom_point(aes(color = Variant_Category), alpha = 0.6, position = position_jitter(width = 0.2, height = 0), size = 1) +
  stat_summary(geom = "point", fun = "mean", color = "black", size = 1) +
  geom_errorbar(stat="summary", fun.data="mean_se", fun.args = list(mult = 1.96), width = 0.1, color = "black") +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("PCR CT Value") +
  ggtitle("Health and Environmental Testing Laboratory (Maine)") +
  theme_bw() +
  scale_x_discrete(labels = addline_format(c(summary_alpha, summary_delta))) + 
  scale_y_continuous(trans = "reverse", breaks = seq(10, 40, by = 10), limits = c(43, 8)) +
  theme(axis.title.x=element_blank(), legend.position="none",
        plot.title = element_text(size = 15)) + 
  theme(plot.title = element_text(size = 11)) +
  stat_compare_means(data = all_CT, mapping = aes(x= Variant_Category, y = CT), 
                     comparisons = list(c("Alpha", "Delta")), method = "t.test", label = "p.signif", vjust = 0.25, hide.ns = FALSE) 

p 
plot_ME_HETL_Jax_Lab1 <- p

mu <- ddply(all_CT, .(Variant_Category), summarize, grp.mean=mean(CT))
p  <- ggplot(all_CT, aes(x=CT, fill = Variant_Category)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = customPalette) +  
  scale_fill_manual(values = customPalette) +
  ylab("Density") +
  xlab("PCR CT Value") +
  ggtitle("Health and Environmental Testing Laboratory (Maine)") +
  theme_bw() +
  labs(color = "Variant Category", fill = "Variant Category") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Variant_Category),
             linetype="dashed")
p
plot_ME_HETL_Jax_Lab1_dens <- p


# Combo Plot --------------------------------------------------------------

p <- ggarrange(plot_CT_Yale_Alpha_Time,
               plot_CT_Yale_Delta_Time,
               labels = c("A", "B"))

ggsave(paste(path_figures, "Fig_S7.pdf", sep=""), p, height = 5, width =10)

p <- ggarrange(plot_CT_Yale,
               plot_CT_Jax,
               plot_MA_MGB,
               plot_ME_HETL_Jax_Lab1)

ggsave(paste(path_figures, "Fig_4.pdf", sep=""), p, height = 7, width =10)

library(gridExtra)
library(ggpubr)

# for grabbing legend for use in a combination plot
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

mylegend <- g_legend(plot_CT_Yale_dens  + theme(legend.position='bottom'))
p <- grid.arrange(ggarrange(plot_CT_Yale_dens + theme(legend.position="none") + xlim(10, 35),
                            plot_CT_Jax_dens + theme(legend.position="none") + xlim(10, 35),
                            plot_MA_MGB_dens + theme(legend.position="none") + xlim(10, 35),
                            plot_ME_HETL_Jax_Lab1_dens + theme(legend.position="none") + xlim(10, 35),
                            nrow = 2, ncol = 2),
                  mylegend,
                  nrow = 2,
                  heights = c(9.5, 0.5))
ggsave(paste(path_figures, "Fig_S10.pdf", sep=""), p, height = 7, width =10)

