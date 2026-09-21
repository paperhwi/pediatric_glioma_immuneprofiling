# Survival analysis: Kaplan-Meier by ecotype + multivariable Cox
# n = 251 main cohort samples with OS data

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

ROOT <- "C:/Users/scott/Downloads/Open PBTA"
OUT  <- file.path(ROOT, "output")
FIG  <- file.path(OUT, "figs_survival")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

# -------- Load ecotype assignment + clinical metadata --------
clus <- fread(file.path(OUT, "ecotype_assignment_k3_annotated.tsv"))

# Subset to samples with OS data
clus$OS_days_num <- as.numeric(clus$OS_days)
clus$OS_event <- as.integer(clus$OS_status == "DECEASED")
surv_df <- clus[!is.na(clus$OS_days_num) & !is.na(clus$OS_event) &
                clus$OS_days_num > 0, ]
cat(sprintf("Samples with OS data: %d / %d\n", nrow(surv_df), nrow(clus)))

# Ensure ecotype is ordered factor for KM legend
surv_df$ecotype <- factor(surv_df$ecotype,
                          levels = c("Inflamed", "Intermediate", "Immune-desert"))

# -------- 1. Overall KM by ecotype --------
cat("\n=== [1] KM by ecotype (overall) ===\n")
fit <- survfit(Surv(OS_days_num, OS_event) ~ ecotype, data = surv_df)
print(fit)
ld <- survdiff(Surv(OS_days_num, OS_event) ~ ecotype, data = surv_df)
lr_p <- 1 - pchisq(ld$chisq, df = length(ld$n) - 1)
cat(sprintf("\nLog-rank p = %.4f\n", lr_p))

# Plot
ec_colors <- c("Inflamed" = "#3B82F6",
               "Intermediate" = "#F59E0B",
               "Immune-desert" = "#9CA3AF")

p_km <- ggsurvplot(
  fit, data = surv_df,
  pval = TRUE, pval.method = TRUE,
  conf.int = FALSE,
  risk.table = TRUE, risk.table.height = 0.28,
  palette = unname(ec_colors[levels(surv_df$ecotype)]),
  legend.labs = levels(surv_df$ecotype),
  legend.title = "Ecotype",
  xlab = "Days from diagnosis", ylab = "Overall survival probability",
  title = sprintf("Overall survival by immune ecotype (n = %d)", nrow(surv_df)),
  ggtheme = theme_bw(base_size = 12)
)
png(file.path(FIG, "KM_overall_by_ecotype.png"), width = 7.5, height = 7, units = "in", res = 200)
print(p_km)
dev.off()
cat("  Wrote KM_overall_by_ecotype.png\n")

# Median survival
med_surv <- summary(fit)$table[, "median", drop = FALSE]
cat("\nMedian OS (days) per ecotype:\n")
print(med_surv)
fwrite(data.frame(ecotype = rownames(med_surv), median_OS_days = med_surv[, 1]),
       file.path(OUT, "KM_median_OS_by_ecotype.tsv"))

# -------- 2. KM stratified by cohort_group --------
cat("\n=== [2] KM by ecotype within each cohort_group ===\n")
cohort_groups <- c("DMG_K27", "DHG_G34", "pHGG_WT", "IHG")

km_strata_p <- data.frame()
for (g in cohort_groups) {
  sub <- surv_df[surv_df$cohort_group == g, ]
  if (nrow(sub) < 10 || length(unique(sub$ecotype)) < 2) {
    cat(sprintf("  [%s] n=%d — skipped (insufficient sample/strata)\n", g, nrow(sub)))
    next
  }
  cat(sprintf("\n  [%s] n=%d\n", g, nrow(sub)))
  print(table(sub$ecotype, sub$OS_event))
  fit_g <- survfit(Surv(OS_days_num, OS_event) ~ ecotype, data = sub)
  ld_g <- survdiff(Surv(OS_days_num, OS_event) ~ ecotype, data = sub)
  p_g <- 1 - pchisq(ld_g$chisq, df = length(ld_g$n) - 1)
  cat(sprintf("    Log-rank p = %.4f\n", p_g))
  km_strata_p <- rbind(km_strata_p,
                       data.frame(cohort_group = g, n = nrow(sub), logrank_p = p_g))
  pp <- tryCatch(ggsurvplot(
    fit_g, data = sub,
    pval = TRUE, conf.int = FALSE,
    risk.table = TRUE, risk.table.height = 0.32,
    palette = unname(ec_colors[levels(droplevels(sub$ecotype))]),
    legend.labs = levels(droplevels(sub$ecotype)),
    legend.title = "Ecotype",
    xlab = "Days from diagnosis", ylab = "OS probability",
    title = sprintf("%s (n=%d)", g, nrow(sub)),
    ggtheme = theme_bw(base_size = 11)
  ), error = function(e) NULL)
  if (!is.null(pp)) {
    png(file.path(FIG, sprintf("KM_by_ecotype_%s.png", g)),
        width = 6.5, height = 6.5, units = "in", res = 200)
    print(pp); dev.off()
    cat(sprintf("    Wrote KM_by_ecotype_%s.png\n", g))
  }
}
fwrite(km_strata_p, file.path(OUT, "KM_stratified_logrank_pvalues.tsv"))

# -------- 3. Multivariable Cox model --------
cat("\n=== [3] Multivariable Cox model ===\n")
# Variables: ecotype + cohort_group + age + sex
surv_df$cohort_group <- factor(surv_df$cohort_group, levels = c("pHGG_WT","DMG_K27","DHG_G34","IHG"))
surv_df$reported_gender <- factor(ifelse(surv_df$reported_gender %in% c("Male","Female"),
                                          surv_df$reported_gender, NA),
                                   levels = c("Male","Female"))

cox_dat <- surv_df[!is.na(reported_gender) & !is.na(age_years), ]
cat(sprintf("Cox model n = %d\n", nrow(cox_dat)))

cox_fit <- tryCatch({
  coxph(Surv(OS_days_num, OS_event) ~ ecotype + cohort_group + age_years + reported_gender,
        data = cox_dat)
}, error = function(e) { print(e); NULL })

if (!is.null(cox_fit)) {
  cat("\n  Cox model summary:\n")
  print(summary(cox_fit))
  # Save coefficients table
  cox_sum <- summary(cox_fit)$coefficients
  cox_ci  <- summary(cox_fit)$conf.int
  cox_tbl <- data.frame(
    term = rownames(cox_sum),
    HR = round(cox_sum[, "exp(coef)"], 3),
    HR_lower95 = round(cox_ci[, "lower .95"], 3),
    HR_upper95 = round(cox_ci[, "upper .95"], 3),
    p_value = signif(cox_sum[, "Pr(>|z|)"], 3)
  )
  fwrite(cox_tbl, file.path(OUT, "Cox_multivariable_coefficients.tsv"))
  cat("\n  Saved Cox_multivariable_coefficients.tsv\n")

  # Forest plot
  png(file.path(FIG, "Cox_forest_multivariable.png"), width = 8, height = 6, units = "in", res = 200)
  print(ggforest(cox_fit, data = cox_dat, fontsize = 0.9))
  dev.off()
  cat("  Wrote Cox_forest_multivariable.png\n")
}

# -------- 4. Univariable HR for ecotype --------
cat("\n=== [4] Univariable Cox for ecotype ===\n")
cox_uni <- coxph(Surv(OS_days_num, OS_event) ~ ecotype, data = surv_df)
print(summary(cox_uni))
uni_sum <- summary(cox_uni)$coefficients
uni_ci  <- summary(cox_uni)$conf.int
uni_tbl <- data.frame(
  term = rownames(uni_sum),
  HR = round(uni_sum[, "exp(coef)"], 3),
  HR_lower95 = round(uni_ci[, "lower .95"], 3),
  HR_upper95 = round(uni_ci[, "upper .95"], 3),
  p_value = signif(uni_sum[, "Pr(>|z|)"], 3)
)
fwrite(uni_tbl, file.path(OUT, "Cox_univariable_ecotype.tsv"))

cat("\nDone.\n")
