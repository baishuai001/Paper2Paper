#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(ggplot2)
  library(survival)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: run_figure4_permutations.R OUTPUT_ROOT PERMUTATIONS THREADS")
out <- normalizePath(args[[1]], mustWork = TRUE)
permutations <- as.integer(args[[2]])
threads <- as.integer(args[[3]])
if (permutations != 5000L) stop("Anchor-aligned Figure 4 requires exactly 5,000 permutations")
tables <- file.path(out, "results", "tables")
figures <- file.path(out, "results", "figures")
audit <- file.path(out, "audit")
inputs <- readRDS(file.path(out, "results", "figure4_survival_inputs.rds"))
observed <- fread(file.path(tables, "Figure4_all_Cox_models.tsv"))
RNGkind("L'Ecuyer-CMRG")
set.seed(20260820)

fit_p <- function(time, event, exposure, clinical, model) {
  dat <- copy(clinical)
  # data.table resolves bare names in j against columns before the calling
  # environment.  Preserve the function arguments under names that cannot be
  # confused with the existing clinical `time` column; otherwise the null
  # analysis silently reuses the observed survival times.
  permuted_time_value <- as.numeric(time)
  event_value <- as.numeric(event)
  exposure_value <- as.numeric(exposure)
  dat[, `:=`(
    permuted_time = permuted_time_value,
    permuted_event = event_value,
    activity_high = exposure_value
  )]
  fit <- tryCatch(
    if (model == "univariable") {
      coxph(Surv(permuted_time, permuted_event) ~ activity_high, data = dat,
            ties = "efron", x = FALSE, y = FALSE, model = FALSE)
    } else {
      coxph(Surv(permuted_time, permuted_event) ~ activity_high + age + sex_model + stage_model + smoking_model,
            data = dat, ties = "efron", x = FALSE, y = FALSE, model = FALSE)
    }, error = function(e) NULL
  )
  if (is.null(fit)) return(1)
  value <- tryCatch(summary(fit)$coefficients["activity_high", "Pr(>|z|)"], error = function(e) 1)
  ifelse(is.finite(value), value, 1)
}

all_results <- list()
summary_rows <- list()
for (input_name in names(inputs)) {
  object <- inputs[[input_name]]
  clinical <- object$clinical
  exposure <- object$activity_high
  for (model_name in c("univariable", "multivariable")) {
    counts <- unlist(mclapply(seq_len(permutations), function(iteration) {
      permuted_time <- sample(clinical$time, replace = FALSE)
      p_values <- vapply(seq_len(nrow(exposure)), function(index) {
        fit_p(permuted_time, clinical$event, as.numeric(exposure[index, ]), clinical, model_name)
      }, numeric(1))
      sum(p_values <= 0.05)
    }, mc.cores = threads, mc.preschedule = TRUE, mc.set.seed = TRUE), use.names = FALSE)
    observed_count <- observed[cohort == object$cohort & endpoint == object$endpoint &
                                 model == model_name, sum(p_value <= 0.05)]
    empirical_p <- mean(counts >= observed_count)
    null_mean <- mean(counts)
    summary_rows[[length(summary_rows) + 1L]] <- data.table(
      cohort = object$cohort, endpoint = object$endpoint, model = model_name,
      observed_significant_TFs = observed_count, null_mean = null_mean,
      null_sd = sd(counts), empirical_p = empirical_p,
      fold_enrichment = ifelse(null_mean > 0, observed_count / null_mean, NA_real_),
      permutations = permutations
    )
    all_results[[length(all_results) + 1L]] <- data.table(
      cohort = object$cohort, endpoint = object$endpoint, model = model_name,
      permutation = seq_len(permutations), significant_TFs = counts,
      observed_significant_TFs = observed_count
    )
  }
}
null_table <- rbindlist(all_results)
summary_table <- rbindlist(summary_rows)
if (any(!is.finite(summary_table$null_sd)) || any(summary_table$null_sd <= 0)) {
  stop("Permutation null is degenerate; verify that shuffled survival times entered the Cox models")
}
fwrite(null_table, file.path(tables, "SupplementaryFigure7_permutation_null_counts.tsv.gz"), sep = "\t")
fwrite(summary_table, file.path(tables, "SupplementaryFigure7_permutation_summary.tsv"), sep = "\t")

null_table[, label := sprintf("%s %s %s", sub("_LUAD$", "", cohort), endpoint, model)]
plot <- ggplot(null_table, aes(significant_TFs)) +
  geom_histogram(binwidth = 1, boundary = -0.5, fill = "#8DA0CB", color = "white") +
  geom_vline(aes(xintercept = observed_significant_TFs), color = "#B2182B", linewidth = 0.7) +
  facet_wrap(~label, scales = "free_y", ncol = 4) +
  labs(title = "Supplementary Figure 7 — 5,000 survival-time permutations",
       subtitle = "Red line: observed number of HC-TFs with Cox p <= 0.05",
       x = "Significant HC-TFs per permutation", y = "Frequency") +
  theme_classic(base_size = 8.2) + theme(strip.text = element_text(size = 7))
ggsave(file.path(figures, "SupplementaryFigure7_permutation_nulls.pdf"), plot,
       width = 12.5, height = 6.2, useDingbats = FALSE)
ggsave(file.path(figures, "SupplementaryFigure7_permutation_nulls.png"), plot,
       width = 12.5, height = 6.2, dpi = 220)

write_json(list(
  status = "passed", permutations = permutations, threads = threads,
  permutation_definition = "survival times shuffled; event status, median activity groups and clinical covariates fixed",
  random_seed = 20260820, RNGkind = "L'Ecuyer-CMRG", summary = summary_table
), file.path(audit, "figure4_permutation_receipt.json"), pretty = TRUE, auto_unbox = TRUE)
cat(toJSON(summary_table), "\n")
