rm(list = ls())

#Setup ----
library(httpgd) #view plots in VS Code
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(flexsdm)

hgd()
path = "/maps/epr26/captain_brazil/af_10km/"
sp_info = read.csv(paste0(path, "species_info.csv"), header = T)
n_sp = nrow(sp_info)

#retrieve model outputs and save evaluation metrics
perf_list = vector("list", n_sp)

for(i in seq_len(n_sp)) {
  if(file.exists(paste0(path, "/models/sdm_model_outputs_", i, ".rds"))) {
    sdm_out = readRDS(paste0(path, "/models/sdm_model_outputs_", i, ".rds"))

    if(!is.null(sdm_out[["models"]])) {
      if(!is.null(sdm_out[["models"]]$gam)) {
        m_out = sdm_out[["models"]]$gam
      } else {
        m_out = NULL
      }
    } else if (!is.null(sdm_out[["model"]])) {
        m_out = sdm_out[["model"]]
    } else {
        m_out = NULL
    }

    if(!is.null(m_out) && !identical(m_out, NA) && file.exists(paste0(path, "/preds/pred_gam_", i, ".tif"))) {
      perf_i = m_out$performance %>%
        filter(threshold == "max_sens_spec") %>%
        dplyr::select(thr_value, TPR_mean, TNR_mean, OR_mean, TSS_mean, AUC_mean) %>%
        mutate(species = i) %>%
        pivot_longer(cols = ends_with("_mean"), names_to = "metric", values_to = "value") %>%
        relocate(species, .before = thr_value) %>%
        mutate(status = "completed")
      cat("Completed saving outputs for species", i, "\n")
    } else {
      cat("Skipping species", i, "as no valid model output\n")
      perf_i = data.frame(species = i, thr_value = NA, metric = NA, value = NA, status = "invalid")
    }

  } else {
    cat("Skipping species", i, "as model not found\n")
    perf_i = data.frame(species = i, thr_value = NA, metric = NA, value = NA, status = "not found")
  }
  perf_list[[i]] = perf_i
}

perf_df = bind_rows(perf_list)

#exclude species based on performance thresholds
sp_excl_tss = subset(perf_df, metric == "TSS_mean" & value < 0.5)$species
sp_excl_tpr = subset(perf_df, metric == "TPR_mean" & value < 0.7)$species
sp_excl_tnr = subset(perf_df, metric == "TNR_mean" & value < 0.7)$species
sp_excl_auc = subset(perf_df, metric == "AUC_mean" & value < 0.8)$species
sp_excl = Reduce(union, list(sp_excl_tss, sp_excl_tpr, sp_excl_tnr, sp_excl_auc))
perf_df$sp_excl = ifelse(perf_df$species %in% sp_excl, T, F)

write.csv(perf_df, paste0(path, "model_performance_all.csv"), row.names = F)

perf_df = read.csv(paste0(path, "model_performance_all.csv"), header = T)

perf_df_plot = subset(perf_df, !is.na(metric)) %>%
  mutate(group = "All species")
perf_df_plot$metric = case_match(perf_df_plot$metric,
                                 "TPR_mean" ~ "Sensitivity",
                                 "TNR_mean" ~ "Specificity",
                                 "OR_mean" ~ "Omission Rate",
                                 "TSS_mean" ~ "True Skill Statistic",
                                 "AUC_mean" ~ "Area Under Curve") %>%
  factor(levels = c("Sensitivity", "Specificity", "Omission Rate", "True Skill Statistic", "Area Under Curve"))
perf_df_plot_sel = subset(perf_df_plot, !sp_excl) %>%
  mutate(group = "Retained species")
perf_df_plot_all = rbind(perf_df_plot, perf_df_plot_sel)

annot_df = perf_df_plot_sel %>%
  group_by(metric) %>%
  summarise(median = median(value, na.rm = T),
            p05 = quantile(value, 0.05, na.rm = T),
            p95 = quantile(value, 0.95, na.rm = T)) %>%
  ungroup() %>%
  mutate(group = "Retained species",
         x = 1,
         y = case_match(metric,
                        "Sensitivity" ~ 0.5,
                        "Specificity" ~ 0.5,
                        "Omission Rate" ~ 0.5,
                        "True Skill Statistic" ~ 0.25,
                        "Area Under Curve" ~ 0.7),
         label = paste0(metric, ": ", round(median, 2), " [", round(p05, 2), "-", round(p95, 2), "]"))

plot_perf = ggplot(perf_df_plot_all, aes(y = value)) +
  geom_violin(aes(x = 1), trim = T, fill = "lightblue", alpha = 0.4) +
  geom_boxplot(aes(x = 1), width = 0.3, outlier.size = 0.5, fill = "white") +
  geom_text(data = annot_df, aes(group = group, x = x, y = y, label = label), size = 5) +
  facet_grid(vars(metric), vars(group), scales = "free_y", switch = "y") +
  labs(title = "", x = "", y = "") +
  theme_bw() +
  theme(strip.placement = "outside",
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 14),
        axis.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12))
plot_perf
ggsave(filename = paste0(path, "plot_model_performance_gam.png"),
       plot = plot_perf, width = 12, height = 12, units = "in", dpi = 300)

#summary table of model performance
perf_df_plot_sel_wide = perf_df_plot_sel %>%
  dplyr::select(-c(thr_value, status, sp_excl, group)) %>%
  pivot_wider(names_from = metric, values_from = value)


completed_id = subset(perf_df, status == "completed")$species %>% unique() #4741 completed species
retained_id = subset(perf_df, status == "completed" & !sp_excl)$species %>% unique() #3885 species not excluded
write.table(retained_id, paste0(path, "species_retained.txt"), sep = "\n", row.names = F, col.names = F)