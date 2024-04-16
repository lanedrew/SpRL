library(readr)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(patchwork)

#### RL ####
cbbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

dens_comp1_rl <- read_csv("./code/results_processing/dens_comp_RL.csv") %>%
  mutate(Replicate = factor(type,
                            levels = c("Observed", "Skew_t", "Skew_Normal", "Normal", "MLR"),
                            labels = c("Observed", "Skewed t", "Skew Normal", "Normal", "MLR")))
dens_comp1_rl |> ggplot(aes(x = growth, group = as.factor(Replicate), color = Replicate)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  labs(x = TeX("Annual Growth $(m^3)$"), y = "Density",
       title = "Record Linkage",
       color = "Replicate Type") +
  theme(text = element_text(family = "serif", size = 14)) -> rl_dens_plot


G_stats_rl <- c(3.287, 1.524, -2.969, 11.649)

summ_stats_skew_t_rl <- read_csv("./code/results_processing/summ_stats_skew_t_rl.csv")

par(mfrow = c(2,2))
hist(summ_stats_skew_t_rl$V1, main = "RL Skew t Model: Mean", xlab = "")
abline(v = G_stats_rl[1], col = "red")

hist(summ_stats_skew_t_rl$V2,  main = "RL Skew t Model: Median", xlab = "")
abline(v = G_stats_rl[2], col = "red")

hist(summ_stats_skew_t_rl$`10%`,  main = "RL Skew t Model: 10th Quantile", xlab = "")
abline(v = G_stats_rl[3], col = "red")

hist(summ_stats_skew_t_rl$`90%`,  main = "RL Skew t Model: 90th Quantile", xlab = "")
abline(v = G_stats_rl[4], col = "red")

summ_stats_skew_n_rl <- read_csv("./code/results_processing/summ_stats_skew_n_rl.csv")

par(mfrow = c(2,2))
hist(summ_stats_skew_n_rl$V1, main = "RL Skew Normal Model: Mean", xlab = "")
abline(v = G_stats_rl[1], col = "red")

hist(summ_stats_skew_n_rl$V2,  main = "RL Skew Normal Model: Median", xlab = "")
abline(v = G_stats_rl[2], col = "red")

hist(summ_stats_skew_n_rl$`10%`,  main = "RL Skew Normal Model: 10th Quantile", xlab = "")
abline(v = G_stats_rl[3], col = "red")

hist(summ_stats_skew_n_rl$`90%`,  main = "RL Skew Normal Model: 90th Quantile", xlab = "")
abline(v = G_stats_rl[4], col = "red")


summ_stats_normal_rl<- read_csv("./code/results_processing/summ_stats_normal_rl.csv")

par(mfrow = c(2,2))
hist(summ_stats_normal_rl$V1, main = "RL Normal Model: Mean", xlab = "")
abline(v = G_stats_rl[1], col = "red")

hist(summ_stats_normal_rl$V2,  main = "RL Normal Model: Median", xlab = "")
abline(v = G_stats_rl[2], col = "red")

hist(summ_stats_normal_rl$`10%`,  main = "RL Normal Model: 10th Quantile", xlab = "")
abline(v = G_stats_rl[3], col = "red")

hist(summ_stats_normal_rl$`90%`,  main = "RL Normal Model: 90th Quantile", xlab = "")
abline(v = G_stats_rl[4], col = "red")


summ_stats_MLR_rl<- read_csv("./code/results_processing/summ_stats_MLR_rl.csv")

# test <- calculate_crps(RL_data$y_rep1, RL_data$y_rep2, G[gc_index])
# RL_data <- readRDS("./code/results_processing/growth_res_MLR_RL.RDS")
# 
# summ_stats_MLR <- apply(RL_data$y_rep1, 1, function(x) c(mean(x), median(x), quantile(x, c(.1, .9))), simplify = FALSE)
# summ_stats_MLR <- do.call(rbind, summ_stats_MLR)
# 
# write_csv(as.data.frame(summ_stats_MLR), "./code/results_processing/summ_stats_MLR_rl.csv")

par(mfrow = c(2,2))
hist(summ_stats_MLR_rl$V1, main = "RL MLR Model: Mean", xlab = "")
abline(v = G_stats_rl[1], col = "red")

hist(summ_stats_MLR_rl$V2,  main = "RL MLR Model: Median", xlab = "")
abline(v = G_stats_rl[2], col = "red")

hist(summ_stats_MLR_rl$`10%`,  main = "RL MLR Model: 10th Quantile", xlab = "")
abline(v = G_stats_rl[3], col = "red")

hist(summ_stats_MLR_rl$`90%`,  main = "RL MLR Model: 90th Quantile", xlab = "")
abline(v = G_stats_rl[4], col = "red")





#### POM ####

dens_comp1_pom <- read_csv("./code/results_processing/dens_comp_pom.csv") %>%
  mutate(Replicate = factor(type,
                            levels = c("Observed", "Skew_t", "Skew_Normal", "Normal", "MLR"),
                            labels = c("Observed", "Skewed t", "Skew Normal", "Normal", "MLR")))
dens_comp1_pom |> ggplot(aes(x = growth, group = as.factor(Replicate), color = Replicate)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  labs(x = TeX("Annual Growth $(m^3)$"), y = "Density",
       title = "Polygon Overlap Matching",
       color = "Replicate Type") +
  theme(text = element_text(family = "serif", size = 14)) -> pom_dens_plot


G_stats_pom <- c(-.81, -.577, -8.68, 6.64)

summ_stats_skew_t_pom <- read_csv("./code/results_processing/summ_stats_skew_t_pom.csv")

par(mfrow = c(2,2))
hist(summ_stats_skew_t_pom$V1, main = "POM Skew t Model: Mean", xlab = "")
abline(v = G_stats_pom[1], col = "red")

hist(summ_stats_skew_t_pom$V2,  main = "POM Skew t Model: Median", xlab = "")
abline(v = G_stats_pom[2], col = "red")

hist(summ_stats_skew_t_pom$`10%`,  main = "POM Skew t Model: 10th Quantile", xlab = "")
abline(v = G_stats_pom[3], col = "red")

hist(summ_stats_skew_t_pom$`90%`,  main = "POM Skew t Model: 90th Quantile", xlab = "")
abline(v = G_stats_pom[4], col = "red")

summ_stats_skew_n_pom <- read_csv("./code/results_processing/summ_stats_skew_n_pom.csv")

par(mfrow = c(2,2))
hist(summ_stats_skew_n_pom$V1, main = "POM Skew Normal Model: Mean", xlab = "")
abline(v = G_stats_pom[1], col = "red")

hist(summ_stats_skew_n_pom$V2,  main = "POM Skew Normal Model: Median", xlab = "")
abline(v = G_stats_pom[2], col = "red")

hist(summ_stats_skew_n_pom$`10%`,  main = "POM Skew Normal Model: 10th Quantile", xlab = "")
abline(v = G_stats_pom[3], col = "red")

hist(summ_stats_skew_n_pom$`90%`,  main = "POM Skew Normal Model: 90th Quantile", xlab = "")
abline(v = G_stats_pom[4], col = "red")


summ_stats_normal_pom<- read_csv("./code/results_processing/summ_stats_normal_pom.csv")

par(mfrow = c(2,2))
hist(summ_stats_normal_pom$V1, main = "POM Normal Model: Mean", xlab = "")
abline(v = G_stats_pom[1], col = "red")

hist(summ_stats_normal_pom$V2,  main = "POM Normal Model: Median", xlab = "")
abline(v = G_stats_pom[2], col = "red")

hist(summ_stats_normal_pom$`10%`,  main = "POM Normal Model: 10th Quantile", xlab = "")
abline(v = G_stats_pom[3], col = "red")

hist(summ_stats_normal_pom$`90%`,  main = "POM Normal Model: 90th Quantile", xlab = "")
abline(v = G_stats_pom[4], col = "red")


summ_stats_MLR_pom<- read_csv("./code/results_processing/summ_stats_MLR_pom.csv")

par(mfrow = c(2,2))
hist(summ_stats_MLR_pom$V1, main = "POM MLR Model: Mean", xlab = "")
abline(v = G_stats_pom[1], col = "red")

hist(summ_stats_MLR_pom$V2,  main = "POM MLR Model: Median", xlab = "")
abline(v = G_stats_pom[2], col = "red")

hist(summ_stats_MLR_pom$`10%`,  main = "POM MLR Model: 10th Quantile", xlab = "")
abline(v = G_stats_pom[3], col = "red")

hist(summ_stats_MLR_pom$`90%`,  main = "POM MLR Model: 90th Quantile", xlab = "")
abline(v = G_stats_pom[4], col = "red")



#### NDM ####


dens_comp1_ndm <- read_csv("./code/results_processing/dens_comp_ndm.csv") %>%
  mutate(Replicate = factor(type,
                            levels = c("Observed", "Skew_t", "Skew_Normal", "Normal", "MLR"),
                            labels = c("Observed", "Skewed t", "Skew Normal", "Normal", "MLR")))

dens_comp1_ndm |> ggplot(aes(x = growth, group = as.factor(Replicate), color = Replicate)) +
  geom_density() +
  theme_bw() +
  scale_color_manual(values = cbbPalette) +
  labs(x = TeX("Annual Growth $(m^3)$"), y = "Density",
       title = "Nearest Distance Matching",
       color = "Replicate Type") +
  theme(text = element_text(family = "serif", size = 14)) -> ndm_dens_plot



G_stats_ndm <- c(4.379, 2.073, -2.72, 14.14)

summ_stats_skew_t_ndm <- read_csv("./code/results_processing/summ_stats_skew_t_ndm.csv")

par(mfrow = c(2,2))
hist(summ_stats_skew_t_ndm$V1, main = "NDM Skew t Model: Mean", xlab = "")
abline(v = G_stats_ndm[1], col = "red")

hist(summ_stats_skew_t_ndm$V2,  main = "NDM Skew t Model: Median", xlab = "")
abline(v = G_stats_ndm[2], col = "red")

hist(summ_stats_skew_t_ndm$`10%`,  main = "NDM Skew t Model: 10th Quantile", xlab = "")
abline(v = G_stats_ndm[3], col = "red")

hist(summ_stats_skew_t_ndm$`90%`,  main = "NDM Skew t Model: 90th Quantile", xlab = "")
abline(v = G_stats_ndm[4], col = "red")


summ_stats_skew_n_ndm <- read_csv("./code/results_processing/summ_stats_skew_n_ndm.csv")

par(mfrow = c(2,2))
hist(summ_stats_skew_n_ndm$V1, main = "NDM Skew Normal Model: Mean", xlab = "")
abline(v = G_stats_ndm[1], col = "red")

hist(summ_stats_skew_n_ndm$V2,  main = "NDM Skew Normal Model: Median", xlab = "")
abline(v = G_stats_ndm[2], col = "red")

hist(summ_stats_skew_n_ndm$`10%`,  main = "NDM Skew Normal Model: 10th Quantile", xlab = "")
abline(v = G_stats_ndm[3], col = "red")

hist(summ_stats_skew_n_ndm$`90%`,  main = "NDM Skew Normal Model: 90th Quantile", xlab = "")
abline(v = G_stats_ndm[4], col = "red")


summ_stats_normal_ndm <- read_csv("./code/results_processing/summ_stats_normal_ndm.csv")

par(mfrow = c(2,2))
hist(summ_stats_normal_ndm$V1, main = "NDM Normal Model: Mean", xlab = "")
abline(v = G_stats_ndm[1], col = "red")

hist(summ_stats_normal_ndm$V1,  main = "NDM Normal Model: Median", xlab = "")
abline(v = G_stats_ndm[2], col = "red")

hist(summ_stats_normal_ndm$`10%`,  main = "NDM Normal Model: 10th Quantile", xlab = "")
abline(v = G_stats_ndm[3], col = "red")

hist(summ_stats_normal_ndm$`90%`,  main = "NDM Normal Model: 90th Quantile", xlab = "")
abline(v = G_stats_ndm[4], col = "red")


summ_stats_MLR_ndm<- read_csv("./code/results_processing/summ_stats_MLR_ndm.csv")

par(mfrow = c(2,2))
hist(summ_stats_MLR_ndm$V1, main = "NDM MLR Model: Mean", xlab = "")
abline(v = G_stats_ndm[1], col = "red")

hist(summ_stats_MLR_ndm$V2,  main = "NDM MLR Model: Median", xlab = "")
abline(v = G_stats_ndm[2], col = "red")

hist(summ_stats_MLR_ndm$`10%`,  main = "NDM MLR Model: 10th Quantile", xlab = "")
abline(v = G_stats_ndm[3], col = "red")

hist(summ_stats_MLR_ndm$`90%`,  main = "NDM MLR Model: 90th Quantile", xlab = "")
abline(v = G_stats_ndm[4], col = "red")



all_dens_plot <- rl_dens_plot + pom_dens_plot + ndm_dens_plot + guide_area() + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'a', tag_suffix = ")")

ggsave(filename = "all_rep_dens_comp_plot.png", plot = all_dens_plot, path = "./plots/F23/",
       width = 20, height = 20, units = "cm", dpi = "retina")
