rm(list = ls())

#Setup ----
library(httpgd) #view plots in VS Code
library(parallel)
library(future) #parallelise lapply() : future_lapply()
library(future.apply) #parallelise lapply(): future_lapply()
library(magrittr)
library(tidyverse)
library(terra)
library(tidyterra)
library(flexsdm)

hgd()
save_path = "/maps/epr26/sdm_captain_out/"

#Set optional user-selected project(s) to run
args = commandArgs(trailingOnly = T)
samp_size_df = read.csv(paste0(save_path, "species_sample_size.csv"), header = T)
n_sp = nrow(samp_size_df)
keep_vars = c(1, 2, 7, 12, 15, 18, 19) #indices of selected bioclim variables after collinearity test

if (length(args) == 2) {
  n0 = as.numeric(args[1])
  n1 = as.numeric(args[2])
} else {
  stop("Error: incorrect number of arguments")
}
cat("Saving outputs for species from", n0, "to", n1, "\n")

#check range coverage of background points
range_coverage_df = read.csv(paste0(save_path, "range_coverage.csv"), header = T)
range_coverage_df_old = range_coverage_df
range_coverage_df$bg_done[1:100] = 0
for(i in seq(n0, n1, 1)) {
  sdm_out = readRDS(paste0(save_path, "/models/sdm_model_outputs_", i, ".rds"))
  if(!is.null(sdm_out$range_coverage)) {
    range_coverage = sdm_out$range_coverage
  } else {
    range_coverage = data.frame(matrix(NA, nrow = 1, ncol = length(keep_vars)))
    rownames(range_coverage) = i
    colnames(range_coverage) = paste0("wc2.1_5m_bio_", keep_vars)
  }
  if(range_coverage_df$bg_done[i] == 0) { #if not done yet
    range_coverage_df[i, ] = c(range_coverage, 1) #save values and mark as done
  }
}

old_min_coverage = apply(range_coverage_df_old[1:100, 1:7], 1, min)
new_min_coverage = apply(range_coverage_df[1:100, 1:7], 1, min)
min_coverage_df = data.frame(min_coverage = c(old_min_coverage, new_min_coverage),
                             samp_size = rep(c("1000", "10000"), each = 100))
plot_min_coverage = ggplot(data = min_coverage_df, aes(x = samp_size, y = min_coverage)) +
  geom_boxplot() +
  labs(title = "Minimum range coverage of background points across selected bioclim variables",
       x = "Sample size", y = "Minimum coverage") +
  theme_minimal()
ggsave(filename = paste0(save_path, "plot_sample_size_range_coverage.png"),
       plot = plot_min_coverage, width = 6, height = 8, units = "in", dpi = 300)
write.csv(range_coverage_df, paste0(save_path, "range_coverage.csv"), row.names = F)


perf_list = vector("list", n1 - n0 + 1)
for(i in seq(n0, n1, 1)) {
  sdm_out = readRDS(paste0(save_path, "/models/sdm_model_outputs_", i, ".rds"))
  if(!is.null(sdm_out$sdm_data)) {
    sdm_models = sdm_out$models
    model_names = names(sdm_models)

    perf_df = lapply(seq_along(sdm_models), function(i) {
      perf_x = sdm_models[[i]]$performance %>%
      filter(threshold == "max_sens_spec") %>%
        dplyr::select(thr_value, TPR_mean, TNR_mean, OR_mean, TSS_mean, AUC_mean) %>%
        mutate(model = model_names[i])
      return(perf_x)
    }) %>%
      bind_rows() %>%
      mutate(species = i)
    perf_list[[i - n0 + 1]] = perf_df %>%
      pivot_longer(cols = ends_with("_mean"), names_to = "metric", values_to = "value")
    write.csv(perf_df, paste0(save_path, "model_performance/performance_", i, ".csv"), row.names = F)
    cat("Completed saving outputs for species", i, "\n")
  } else {
    cat("Skipping species", i, "as models not found\n")
    next
  }
}

perf_all = perf_list %>%
  bind_rows()

plot_perf = ggplot(data = perf_all, aes(x = model, y = value)) +
  geom_boxplot() +
  facet_wrap(~ metric, scales = "free_y") +
  labs(title = "Model performance across species 1 - 100",
       x = "Model", y = "Metric value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, size = 16),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))
ggsave(filename = paste0(save_path, "plot_model_performance.png"),
       plot = plot_perf, width = 10, height = 8, units = "in", dpi = 300)


write.csv(perf_all, paste0(save_path, "model_performance/performance_all.csv"), row.names = F)

