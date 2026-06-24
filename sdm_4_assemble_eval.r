rm(list = ls())

# Setup ----
library(renv)
library(dplyr)
library(stringr)

on_cluster = nchar(Sys.getenv("SCRATCH_PATH")) > 0

if (on_cluster) {
  cat("=== Running on CLUSTER ===\n"); flush.console()
  
  in_path = paste0(Sys.getenv("PROJECT_PATH"), "/data_in/")
  model_path = paste0(Sys.getenv("PROJECT_PATH"), "/AF_foundmod_sdm_out/")
  out_path = paste0(Sys.getenv("SCRATCH_PATH"), "/pred_endangered_60k/")
  
  cat("Cluster environment detected:\n"); flush.console()
  cat("  PROJECT_PATH:", Sys.getenv("PROJECT_PATH"), "\n"); flush.console()
  cat("  SCRATCH_PATH:", Sys.getenv("SCRATCH_PATH"), "\n"); flush.console()
  
} else {
  cat("=== Running LOCALLY ===\n"); flush.console()
  cat("  Working directory:", getwd(), "\n"); flush.console()

  in_path = "../DATA/CAPTAIN_Brazil_in/"
  out_path = "../DATA/CAPTAIN_Brazil_out/"
  scratch_path = tempdir()
}

id_done = list.files(out_path) |>
  str_extract("eval_\\d+\\.csv") |>
  na.omit() |>
  str_remove_all("eval_|\\.csv") |>
  as.numeric() |>
  sort()

eval_df = lapply(id_done, function(i) {
  eval = read.csv(paste0(out_path, "eval_", i, ".csv"), header = T)
  return(eval)
}) |>
  bind_rows()

eval_part_list = lapply(id_done, function(i) {
  eval_part = read.csv(paste0(out_path, "eval_part_", i, ".csv"), header = T)
  return(eval_part)
})


eval_df = m_out$performance |>
  filter(threshold == "max_sens_spec") |>
  dplyr::select(thr_value, TPR_mean, TNR_mean, TSS_mean, AUC_mean) |>
  mutate(species = i) |>
  pivot_longer(cols = ends_with("_mean"), names_to = "metric", values_to = "value") |>
  relocate(species, .before = thr_value) |>
  mutate(criterion = sapply(metric, AssignCriterion)) |>
  mutate(exclude = value < criterion)

eval_part_df = m_out$performance_part |>
  filter(threshold == "max_sens_spec") |>
  dplyr::select(thr_value, TPR, TNR, TSS, AUC) |>
  summarize(across(thr_value:AUC, function(x) mean(x, na.rm = T))) |>
  mutate(species = i) |>
  relocate(species, .before = thr_value)
write.csv(eval_df, paste0(out_path, "eval_", i, ".csv"), row.names = F)
write.csv(eval_part_df, paste0(out_path, "eval_part_", i, ".csv"), row.names = F)