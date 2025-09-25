rm(list = ls())

#Setup ----
library(httpgd) #view plots in VS Code
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(flexsdm)

hgd()
save_path = "/maps/epr26/sdm_captain_out/"
samp_size_df = read.csv(paste0(save_path, "species_sample_size.csv"), header = T)
n_sp = nrow(samp_size_df)

#retrieve model outputs and save evaluation metrics
perf_list = vector("list", n_sp)

for(i in seq_len(n_sp)) {
  if(file.exists(paste0(save_path, "/models/sdm_model_outputs_", i, ".rds"))) {
    sdm_out = readRDS(paste0(save_path, "/models/sdm_model_outputs_", i, ".rds"))

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

    if(!is.null(m_out) && !identical(m_out, NA) && file.exists(paste0(save_path, "/preds/pred_gam_", i, ".tif"))) {
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

sp_excl_tss = subset(perf_df, metric == "TSS_mean" & value < 0.5)$species
sp_excl_tpr = subset(perf_df, metric == "TPR_mean" & value < 0.7)$species
sp_excl_tnr = subset(perf_df, metric == "TNR_mean" & value < 0.7)$species
sp_excl_auc = subset(perf_df, metric == "AUC_mean" & value < 0.8)$species

sp_excl = Reduce(union, list(sp_excl_tss, sp_excl_tpr, sp_excl_tnr, sp_excl_auc))
perf_df$sp_excl = ifelse(perf_df$species %in% sp_excl, T, F)

plot_perf = ggplot(data = perf_df, aes(y = value)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "GAM performance across all species",
       x = "", y = "Metric value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(filename = paste0(save_path, "plot_model_performance_all.png"),
       plot = plot_perf, width = 10, height = 8, units = "in", dpi = 300)

plot_perf_excl = ggplot(data = subset(perf_df, species %in% sp_excl == F), aes(y = value)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "GAM performance across species with TSS >= 0.5",
       x = "", y = "Metric value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(filename = paste0(save_path, "plot_model_performance_excl.png"),
       plot = plot_perf_excl, width = 10, height = 8, units = "in", dpi = 300)

write.csv(perf_df, paste0(save_path, "model_performance_all.csv"), row.names = F)

completed_id = subset(perf_df, status == "completed")$species %>% unique() #4741 completed species
retained_id = subset(perf_df, status == "completed" & !sp_excl)$species %>% unique() #3885 species not excluded
write.table(retained_id, paste0(save_path, "species_retained.txt"), sep = "\n", row.names = F, col.names = F)




