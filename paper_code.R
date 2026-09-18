library(dplyr)
library(tidyr)
library(purrr)
library(officer)
library(flextable)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(readxl)
library(geepack)
library(splines)

# ============================================================
# SETUP
# ============================================================
out_dir <- "~/Desktop/CMC/MAL-ED/Anemia analysis/Final/Paper/R_outputs/"
if (!dir.exists(path.expand(out_dir))) dir.create(path.expand(out_dir), recursive = TRUE)

data <- read_excel("~/Desktop/CMC/MAL-ED/Anemia analysis/Final/Paper/Final paper docs/maled_data.xlsx")
data <- as.data.frame(data)

# ============================================================
# STANDARDIZE COLUMN NAMES
# ============================================================
data <- data %>%
  rename(
    pid    = si_no,
    CAFSEX = sex,
    
    Hb7   = Hb_7,   Hb15  = Hb_15,  Hb24  = Hb_24,
    Hb36  = Hb_36,  Hb60  = Hb_60,  Hb84  = Hb_84,
    Hb108 = Hb_108, Hb144 = Hb_144,
    
    Ferritin7 = Ferritin_7, Ferritin15 = Ferritin_15, Ferritin24 = Ferritin_24,
    Ferritin108 = Ferritin_108, Ferritin144 = Ferritin_144,
    
    tfr7 = Tfr_7, tfr15 = Tfr_15, tfr24 = Tfr_24,
    tfr108 = Tfr_108, tfr144 = Tfr_144,
    
    bodyiron7 = Bodyiron_7, bodyiron15 = Bodyiron_15, bodyiron24 = Bodyiron_24,
    bodyiron108 = Bodyiron_108, bodyiron144 = Bodyiron_144,
    
    B12_9y  = B12_108,
    B12_12y = B12_144,
    
    Lead_12y = Lead_144
  )

cols_to_convert <- c(
  "Hb7","Hb15","Hb24","Hb36","Hb60","Hb84","Hb108","Hb144",
  "Ferritin7","Ferritin15","Ferritin24","Ferritin108","Ferritin144",
  "tfr7","tfr15","tfr24","tfr108","tfr144",
  "bodyiron7","bodyiron15","bodyiron24","bodyiron108","bodyiron144",
  "B12_9y","B12_12y",
  "Lead_15","Lead_24","Lead_12y"
)
for (col in cols_to_convert) {
  if (col %in% names(data)) data[[col]] <- as.numeric(as.character(data[[col]]))
}

# ============================================================
# NOTE ON DENOMINATOR CONVENTION (applies throughout this script)
# ============================================================
# All prevalence / completeness tables that report iron, ferritin, TfR,
# body iron, Vitamin B12, and Lead values are restricted to children who
# ALSO have a concurrent hemoglobin (Hb) measurement at that same visit.
# This is enforced in two places:
#   1) classify_row() / long_data below only adds a row for a child-age
#      combination if Hb is non-missing at that age (unchanged from
#      before), which is what Tables 1, 3, 4, and 7 are built from.
#   2) Table 8 (completeness) now explicitly requires the same condition
#      when counting each biomarker's N, so its denominators and
#      numerators are constructed the same way as Tables 1/3/7 use.
# This eliminates the previous inconsistency where Table 8 counted
# TfR/Lead/B12 among ALL children with that assay (regardless of Hb
# status), which could exceed the Hb-based denominator and produce
# impossible (>100%) percentages.
# ============================================================

# ============================================================
# CLASSIFICATION FUNCTIONS (WHO age/sex-specific Hb cut-offs)
# ============================================================
hb_cutoff <- function(hb, age_months, sex = NA) {
  if (is.na(hb)) return(NA)
  if      (age_months <= 23)                       return(hb < 10.5)
  else if (age_months >= 24  && age_months <= 59)  return(hb < 11.0)
  else if (age_months >= 60  && age_months <= 131) return(hb < 11.5)
  else if (age_months >= 132 && age_months <= 167) return(hb < 12.0)
  else { if (!is.na(sex) && sex == "Boys") return(hb < 13.0) else return(hb < 12.0) }
}

hb_severity <- function(hb, age_months, sex = NA) {
  if (is.na(hb)) return(NA_character_)
  if (age_months <= 23) {
    if (hb >= 10.5) "No Anaemia" else if (hb >= 9.5) "Mild" else if (hb >= 7.0) "Moderate" else "Severe"
  } else if (age_months <= 59) {
    if (hb >= 11.0) "No Anaemia" else if (hb >= 10.0) "Mild" else if (hb >= 7.0) "Moderate" else "Severe"
  } else if (age_months <= 131) {
    if (hb >= 11.5) "No Anaemia" else if (hb >= 11.0) "Mild" else if (hb >= 8.0) "Moderate" else "Severe"
  } else if (age_months <= 167) {
    if (hb >= 12.0) "No Anaemia" else if (hb >= 11.0) "Mild" else if (hb >= 8.0) "Moderate" else "Severe"
  } else {
    no_an <- if (!is.na(sex) && sex == "Boys") 13.0 else 12.0
    if (hb >= no_an) "No Anaemia" else if (hb >= 11.0) "Mild" else if (hb >= 8.0) "Moderate" else "Severe"
  }
}

iron_cutoff <- function(iron) ifelse(is.na(iron), NA, iron < 0)

b12_status <- function(b12) {
  if (is.na(b12)) return(NA)
  else if (b12 < 200) return("Deficient")
  else if (b12 < 300) return("Insufficient")
  else return("Sufficient")
}

# ============================================================
# BUILD long_data — used for Tables 1, 3, 4 (ages 7,15,24,108,144)
# NOTE: classify_row() only emits a row when Hb is non-missing at that
# age, so every downstream table built from long_data (Tables 1, 3, 4,
# 7) is, by construction, restricted to children with a concurrent Hb
# measurement. This is intentional and is now the standard convention
# applied consistently across all tables (see note above and Table 8).
# ============================================================
classification_data <- data %>%
  select(pid, CAFSEX,
         Hb7, Hb15, Hb24, Hb108, Hb144,
         bodyiron7, bodyiron15, bodyiron24, bodyiron108, bodyiron144,
         B12_9y, B12_12y,
         Lead_15, Lead_24, Lead_12y) %>%
  mutate(sex_group = case_when(CAFSEX == 1 ~ "Boys", CAFSEX == 2 ~ "Girls", TRUE ~ NA_character_))

classify_row <- function(pid, sex_group,
                         hb7, hb15, hb24, hb108, hb144,
                         bodyiron7, bodyiron15, bodyiron24, bodyiron108, bodyiron144,
                         B12_9y, B12_12y,
                         Lead_15, Lead_24, Lead_12y) {
  
  ages      <- c(7, 15, 24, 108, 144)
  hb_vals   <- c(hb7, hb15, hb24, hb108, hb144)
  iron_vals <- c(bodyiron7, bodyiron15, bodyiron24, bodyiron108, bodyiron144)
  b12_vals  <- c(NA, NA, NA, B12_9y, B12_12y)
  lead_vals <- c(NA, Lead_15, Lead_24, NA, Lead_12y)
  
  results <- tibble()
  for (i in seq_along(ages)) {
    # Restrict to children with a concurrent Hb measurement at this age
    if (!is.na(hb_vals[i])) {
      results <- bind_rows(results, tibble(
        pid            = pid,
        sex_group      = sex_group,
        age_months     = ages[i],
        hb             = hb_vals[i],
        bodyiron       = iron_vals[i],
        b12            = b12_vals[i],
        lead           = lead_vals[i],
        anemic         = hb_cutoff(hb_vals[i], ages[i], sex_group),
        hb_severity    = hb_severity(hb_vals[i], ages[i], sex_group),
        iron_deficient = iron_cutoff(iron_vals[i]),
        b12_category   = b12_status(b12_vals[i])
      ))
    }
  }
  results
}

long_data <- classification_data %>%
  rowwise() %>%
  do(classify_row(
    .$pid, .$sex_group,
    .$Hb7, .$Hb15, .$Hb24, .$Hb108, .$Hb144,
    .$bodyiron7, .$bodyiron15, .$bodyiron24, .$bodyiron108, .$bodyiron144,
    .$B12_9y, .$B12_12y,
    .$Lead_15, .$Lead_24, .$Lead_12y
  )) %>%
  ungroup() %>%
  mutate(
    hb_status   = ifelse(anemic, "Anemic", "Non-Anemic"),
    lead_status = case_when(is.na(lead) ~ NA_character_, lead >= 5 ~ "Lead Elevated", TRUE ~ "Lead Normal")
  )

# ============================================================
# TABLE 1: Iron / Vitamin B12 / Lead prevalence by age & sex
# ============================================================
fmt_pct <- function(n, total) {
  if (is.na(total) || total == 0) return("—")
  paste0(round(100 * n / total, 1), "% (", n, "/", total, ")")
}

compute_row1 <- function(df) {
  total_iron <- sum(!is.na(df$iron_deficient))
  total_b12  <- sum(!is.na(df$b12_category))
  total_lead <- sum(!is.na(df$lead_status))
  tibble(
    Iron_Deficient   = fmt_pct(sum(df$iron_deficient == TRUE,  na.rm = TRUE), total_iron),
    Iron_Sufficient  = fmt_pct(sum(df$iron_deficient == FALSE, na.rm = TRUE), total_iron),
    B12_Deficient    = fmt_pct(sum(df$b12_category == "Deficient",    na.rm = TRUE), total_b12),
    B12_Insufficient = fmt_pct(sum(df$b12_category == "Insufficient", na.rm = TRUE), total_b12),
    B12_Sufficient   = fmt_pct(sum(df$b12_category == "Sufficient",   na.rm = TRUE), total_b12),
    Lead_Normal      = fmt_pct(sum(df$lead_status == "Lead Normal",   na.rm = TRUE), total_lead),
    Lead_Elevated    = fmt_pct(sum(df$lead_status == "Lead Elevated", na.rm = TRUE), total_lead)
  )
}

age_levels1 <- c(7, 15, 24, 108, 144)
age_labels1 <- c("7 months", "15 months", "24 months", "108 months", "144 months")

rows1 <- list()
for (i in seq_along(age_levels1)) {
  age    <- age_levels1[i]
  lbl    <- age_labels1[i]
  age_df <- long_data %>% filter(age_months == age)
  for (sex in c("All Children", "Boys", "Girls")) {
    sub_df <- if (sex == "All Children") age_df else age_df %>% filter(sex_group == sex)
    rows1[[length(rows1) + 1]] <- compute_row1(sub_df) %>%
      mutate(Age = lbl, Sex = sex, .before = everything())
  }
}
table1 <- bind_rows(rows1) %>%
  rename(
    `Deficient` = Iron_Deficient, `Sufficient` = Iron_Sufficient,
    `Deficient ` = B12_Deficient, `Insufficient` = B12_Insufficient, `Sufficient ` = B12_Sufficient,
    `Lead Normal` = Lead_Normal, `Lead Elevated` = Lead_Elevated
  )

table1_caption <- paste0(
  "Table 1: Percentage of iron deficiency, Vitamin B12 deficiency, and Lead exposure among children ",
  "at 7, 15, 24, 108, and 144 months of age. Restricted to children with a concurrent hemoglobin (Hb) ",
  "measurement at each age; denominators therefore match those used in Tables 3 and 7 and in the ",
  "Table 8 completeness figures (see Table 8 footnote)."
)

ft1 <- flextable(table1) %>%
  add_header_row(values = c("", "", "Iron (mg/kg)", "Vitamin B12 (pg/ml)", "Lead (µg/dL)"),
                 colwidths = c(1, 1, 2, 3, 2)) %>%
  merge_v(j = "Age", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  align(align = "center", part = "all") %>%
  align(j = 2, align = "left", part = "body") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  valign(j = 1, valign = "center", part = "body") %>%
  autofit() %>%
  set_caption(table1_caption)

doc1 <- read_docx() %>%
  body_add_par(table1_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft1)
print(doc1, target = paste0(out_dir, "Table1_Iron_B12_Lead_Prevalence.docx"))

# ============================================================
# WILSON 95% CI HELPERS
# Used by Table 10 (already) and by the new Tables 12-14 below, which
# add exact 95% CIs to the small-N proportions in Tables 2, 3, and 4
# per reviewer note #10 (zero-cell / small-N proportions need CIs,
# especially the 108-144 month cells with N < 20 such as the 0/17 vs
# 17/17 lead split at 144 months, Issue #3).
# ============================================================
wilson_ci <- function(x, n, z = 1.96) {
  if (is.na(n) || n == 0) return(c(lo = NA_real_, hi = NA_real_))
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2))) / denom
  c(lo = max(0, center - margin), hi = min(1, center + margin))
}

# Formats "n/total" as a percentage with an appended Wilson 95% CI,
# e.g. "100.0% (17/17) [80.5-100.0]"
fmt_pct_ci <- function(n, total) {
  if (is.na(total) || total == 0) return("—")
  p  <- 100 * n / total
  ci <- 100 * wilson_ci(n, total)
  sprintf("%.1f%% (%d/%d) [%.1f\u2013%.1f]", p, n, total, ci["lo"], ci["hi"])
}

# ============================================================
# TABLE 2: Anaemia severity, 7–144 months, ALL Hb timepoints, all children
# ============================================================
hb_long_all <- data %>%
  select(pid, CAFSEX, Hb7, Hb15, Hb24, Hb36, Hb60, Hb84, Hb108, Hb144) %>%
  pivot_longer(cols = -c(pid, CAFSEX), names_to = "age_label", values_to = "hb") %>%
  mutate(
    age_months = case_when(
      age_label == "Hb7"   ~ 7,   age_label == "Hb15"  ~ 15,  age_label == "Hb24"  ~ 24,
      age_label == "Hb36"  ~ 36,  age_label == "Hb60"  ~ 60,  age_label == "Hb84"  ~ 84,
      age_label == "Hb108" ~ 108, age_label == "Hb144" ~ 144,
      TRUE ~ NA_real_
    ),
    CAFSEX   = factor(CAFSEX, levels = c(1, 2), labels = c("Boys", "Girls")),
    anemic   = mapply(hb_cutoff,   hb, age_months, as.character(CAFSEX)),
    severity = mapply(hb_severity, hb, age_months, as.character(CAFSEX)),
    severity = factor(severity, levels = c("No Anaemia", "Mild", "Moderate", "Severe"))
  ) %>%
  filter(!is.na(hb))

fmt_pct_vec <- function(n, total) {
  ifelse(is.na(total) | total == 0, "—", paste0(round(100 * n / total, 1), "% (", n, "/", total, ")"))
}

severity_summary <- hb_long_all %>%
  group_by(age_months) %>%
  summarise(
    Total_N    = n(),
    n_anemic   = sum(anemic == TRUE, na.rm = TRUE),
    n_mild     = sum(severity == "Mild",     na.rm = TRUE),
    n_moderate = sum(severity == "Moderate", na.rm = TRUE),
    n_severe   = sum(severity == "Severe",   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Age (Months)` = age_months,
    Anemic         = fmt_pct_vec(n_anemic,   Total_N),
    Mild           = fmt_pct_vec(n_mild,     Total_N),
    Moderate       = fmt_pct_vec(n_moderate, Total_N),
    Severe         = fmt_pct_vec(n_severe,   Total_N)
  ) %>%
  arrange(age_months) %>%
  select(`Age (Months)`, Anemic, Mild, Moderate, Severe)

ft2 <- flextable(severity_summary) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  align(align = "center", part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption("Table 2: Anemia severity among children of age 7-144 months using new WHO age-specific hemoglobin cut-offs.")

doc2 <- read_docx() %>%
  body_add_par("Table 2: Anemia severity among children of age 7-144 months using new WHO age-specific hemoglobin cut-offs.", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft2)
print(doc2, target = paste0(out_dir, "Table2_Anaemia_Severity.docx"))

# ============================================================
# TABLE 12: Anaemia severity WITH Wilson 95% CI (companion to Table 2)
# Kept as a separate table rather than modifying Table 2 in place, so
# the original percentage-only presentation is preserved and this CI
# version can be used wherever reviewers want exact uncertainty bounds
# on the small-N severity cells.
# ============================================================
severity_summary_ci <- hb_long_all %>%
  group_by(age_months) %>%
  summarise(
    Total_N    = n(),
    n_anemic   = sum(anemic == TRUE, na.rm = TRUE),
    n_mild     = sum(severity == "Mild",     na.rm = TRUE),
    n_moderate = sum(severity == "Moderate", na.rm = TRUE),
    n_severe   = sum(severity == "Severe",   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    `Age (Months)` = age_months,
    Anemic         = mapply(fmt_pct_ci, n_anemic,   Total_N),
    Mild           = mapply(fmt_pct_ci, n_mild,     Total_N),
    Moderate       = mapply(fmt_pct_ci, n_moderate, Total_N),
    Severe         = mapply(fmt_pct_ci, n_severe,   Total_N)
  ) %>%
  arrange(age_months) %>%
  select(`Age (Months)`, Anemic, Mild, Moderate, Severe)

ft12 <- flextable(severity_summary_ci) %>%
  bold(part = "header") %>%
  fontsize(size = 8, part = "all") %>%
  align(align = "center", part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption("Table 12: Anemia severity among children of age 7-144 months using new WHO age-specific hemoglobin cut-offs, with exact Wilson 95% confidence intervals shown in brackets. Companion table to Table 2, recommended for small-N cells (particularly 108 and 144 months).")

doc12 <- read_docx() %>%
  body_add_par("Table 12: Anemia severity by age, with Wilson 95% CI (companion to Table 2).", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft12)
print(doc12, target = paste0(out_dir, "Table12_Anaemia_Severity_WithCI.docx"))

# ============================================================
# TABLE 3: Iron / B12 / Lead prevalence by Hb status (Anemic vs Non-Anemic)
# ============================================================
create_pct_row <- function(data, age_val, hb_status_val) {
  sub <- data %>% filter(age_months == age_val, hb_status == hb_status_val)
  
  iron_sub   <- sub %>% filter(!is.na(iron_deficient)); iron_n <- nrow(iron_sub)
  iron_def_n <- sum(iron_sub$iron_deficient == TRUE,  na.rm = TRUE)
  iron_suf_n <- sum(iron_sub$iron_deficient == FALSE, na.rm = TRUE)
  iron_def   <- if (iron_n > 0) sprintf("%.1f%% (%d/%d)", 100 * iron_def_n / iron_n, iron_def_n, iron_n) else "N/A"
  iron_suf   <- if (iron_n > 0) sprintf("%.1f%% (%d/%d)", 100 * iron_suf_n / iron_n, iron_suf_n, iron_n) else "N/A"
  
  b12_sub   <- sub %>% filter(!is.na(b12_category)); b12_n <- nrow(b12_sub)
  b12_def_n <- sum(b12_sub$b12_category == "Deficient",    na.rm = TRUE)
  b12_ins_n <- sum(b12_sub$b12_category == "Insufficient", na.rm = TRUE)
  b12_suf_n <- sum(b12_sub$b12_category == "Sufficient",   na.rm = TRUE)
  b12_def   <- if (b12_n > 0) sprintf("%.1f%% (%d/%d)", 100 * b12_def_n / b12_n, b12_def_n, b12_n) else "N/A"
  b12_ins   <- if (b12_n > 0) sprintf("%.1f%% (%d/%d)", 100 * b12_ins_n / b12_n, b12_ins_n, b12_n) else "N/A"
  b12_suf   <- if (b12_n > 0) sprintf("%.1f%% (%d/%d)", 100 * b12_suf_n / b12_n, b12_suf_n, b12_n) else "N/A"
  
  has_lead    <- age_val %in% c(15, 24, 144)
  lead_sub    <- sub %>% filter(!is.na(lead)); lead_n <- nrow(lead_sub)
  lead_norm_n <- sum(lead_sub$lead <  5, na.rm = TRUE)
  lead_elev_n <- sum(lead_sub$lead >= 5, na.rm = TRUE)
  lead_norm   <- if (has_lead && lead_n > 0) sprintf("%.1f%% (%d/%d)", 100 * lead_norm_n / lead_n, lead_norm_n, lead_n) else "N/A"
  lead_elev   <- if (has_lead && lead_n > 0) sprintf("%.1f%% (%d/%d)", 100 * lead_elev_n / lead_n, lead_elev_n, lead_n) else "N/A"
  
  tibble(
    Age = age_val, `Hb Status` = hb_status_val,
    `Iron Deficient` = iron_def, `Iron Sufficient` = iron_suf,
    `B12 Deficient` = b12_def, `B12 Insufficient` = b12_ins, `B12 Sufficient` = b12_suf,
    `Lead Normal` = lead_norm, `Lead Elevated` = lead_elev
  )
}

age_vals3 <- c(7, 15, 24, 108, 144)
hb_vals3  <- c("Anemic", "Non-Anemic")

pct_table <- bind_rows(lapply(age_vals3, function(a) bind_rows(lapply(hb_vals3, function(h) create_pct_row(long_data, a, h)))))

age_labels3 <- c(`7` = "7 Month", `15` = "15 Month", `24` = "24 Month", `108` = "108 Month", `144` = "144 Month")
pct_table <- pct_table %>% mutate(Age = age_labels3[as.character(Age)])

ft3 <- flextable(pct_table) %>%
  set_header_labels(
    Age = "Age", `Hb Status` = "Hb Status",
    `Iron Deficient` = "Deficient", `Iron Sufficient` = "Sufficient",
    `B12 Deficient` = "Deficient", `B12 Insufficient` = "Insufficient", `B12 Sufficient` = "Sufficient",
    `Lead Normal` = "Normal", `Lead Elevated` = "Elevated"
  ) %>%
  add_header_row(
    top = TRUE,
    values = c("", "", "Iron (mg/kg)", "Vitamin B12 (pg/ml)", "Lead (µg/dL)"),
    colwidths = c(1, 1, 2, 3, 2)
  ) %>%
  merge_v(j = "Age") %>%
  valign(j = "Age", valign = "top", part = "body") %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "body") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  bold(part = "header") %>%
  autofit()

table3_caption <- paste0(
  "Table 3: Prevalence of iron deficiency, vitamin B12 deficiency and blood lead levels among anemic and ",
  "non-anemic (New cut-off) children at 7, 15, 24, 108, and 144 months of age. Restricted to children with a ",
  "concurrent hemoglobin (Hb) measurement at each age (same convention as Table 1 and Table 8)."
)

doc3 <- read_docx() %>%
  body_add_par(table3_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft3)
print(doc3, target = paste0(out_dir, "Table3_Iron_B12_Lead_by_HbStatus.docx"))

# ============================================================
# TABLE 13: Iron / B12 / Lead prevalence by Hb status WITH Wilson 95% CI
# (companion to Table 3). This is where the 144-month lead split
# (0/17 Non-Anemic... 17/17 Anemic, Issue #3) gets its exact CI, e.g.
# 17/17 = 100.0% [80.5-100.0], and 0/x cells get an honest upper bound
# instead of a bare "0.0%".
# ============================================================
create_pct_row_ci <- function(data, age_val, hb_status_val) {
  sub <- data %>% filter(age_months == age_val, hb_status == hb_status_val)
  
  iron_sub <- sub %>% filter(!is.na(iron_deficient)); iron_n <- nrow(iron_sub)
  iron_def <- if (iron_n > 0) fmt_pct_ci(sum(iron_sub$iron_deficient == TRUE,  na.rm = TRUE), iron_n) else "N/A"
  iron_suf <- if (iron_n > 0) fmt_pct_ci(sum(iron_sub$iron_deficient == FALSE, na.rm = TRUE), iron_n) else "N/A"
  
  b12_sub <- sub %>% filter(!is.na(b12_category)); b12_n <- nrow(b12_sub)
  b12_def <- if (b12_n > 0) fmt_pct_ci(sum(b12_sub$b12_category == "Deficient",    na.rm = TRUE), b12_n) else "N/A"
  b12_ins <- if (b12_n > 0) fmt_pct_ci(sum(b12_sub$b12_category == "Insufficient", na.rm = TRUE), b12_n) else "N/A"
  b12_suf <- if (b12_n > 0) fmt_pct_ci(sum(b12_sub$b12_category == "Sufficient",   na.rm = TRUE), b12_n) else "N/A"
  
  has_lead  <- age_val %in% c(15, 24, 144)
  lead_sub  <- sub %>% filter(!is.na(lead)); lead_n <- nrow(lead_sub)
  lead_norm <- if (has_lead && lead_n > 0) fmt_pct_ci(sum(lead_sub$lead <  5, na.rm = TRUE), lead_n) else "N/A"
  lead_elev <- if (has_lead && lead_n > 0) fmt_pct_ci(sum(lead_sub$lead >= 5, na.rm = TRUE), lead_n) else "N/A"
  
  tibble(
    Age = age_val, `Hb Status` = hb_status_val,
    `Iron Deficient` = iron_def, `Iron Sufficient` = iron_suf,
    `B12 Deficient` = b12_def, `B12 Insufficient` = b12_ins, `B12 Sufficient` = b12_suf,
    `Lead Normal` = lead_norm, `Lead Elevated` = lead_elev
  )
}

pct_table_ci <- bind_rows(lapply(age_vals3, function(a) bind_rows(lapply(hb_vals3, function(h) create_pct_row_ci(long_data, a, h)))))
pct_table_ci <- pct_table_ci %>% mutate(Age = age_labels3[as.character(Age)])

ft13 <- flextable(pct_table_ci) %>%
  set_header_labels(
    Age = "Age", `Hb Status` = "Hb Status",
    `Iron Deficient` = "Deficient", `Iron Sufficient` = "Sufficient",
    `B12 Deficient` = "Deficient", `B12 Insufficient` = "Insufficient", `B12 Sufficient` = "Sufficient",
    `Lead Normal` = "Normal", `Lead Elevated` = "Elevated"
  ) %>%
  add_header_row(
    top = TRUE,
    values = c("", "", "Iron (mg/kg)", "Vitamin B12 (pg/ml)", "Lead (µg/dL)"),
    colwidths = c(1, 1, 2, 3, 2)
  ) %>%
  merge_v(j = "Age") %>%
  valign(j = "Age", valign = "top", part = "body") %>%
  align(align = "center", part = "all") %>%
  align(j = 1:2, align = "left", part = "body") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  bold(part = "header") %>%
  fontsize(size = 7, part = "all") %>%
  autofit()

table13_caption <- paste0(
  "Table 13: Prevalence of iron deficiency, vitamin B12 deficiency and blood lead levels among anemic and ",
  "non-anemic (New cut-off) children at 7, 15, 24, 108, and 144 months of age, with exact Wilson 95% confidence ",
  "intervals shown in brackets. Companion table to Table 3; restricted to children with a concurrent hemoglobin ",
  "(Hb) measurement at each age. Recommended for interpreting small-N cells, e.g. the 144-month lead split among ",
  "anemic children (see Table 3), where the point estimate alone can overstate certainty."
)

doc13 <- read_docx() %>%
  body_add_par(table13_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft13)
print(doc13, target = paste0(out_dir, "Table13_Iron_B12_Lead_by_HbStatus_WithCI.docx"))

# ============================================================
# TABLE 4: Combined iron/B12/lead deficiency by anaemia status, 144 months
# ============================================================
compute_group4 <- function(fd) {
  n      <- nrow(fd)
  lead_n <- sum(!is.na(fd$lead))
  elev_n <- sum(fd$lead >= 5, na.rm = TRUE)
  norm_n <- sum(fd$lead < 5,  na.rm = TRUE)
  elev_txt <- if (lead_n > 0) sprintf("%d (%.2f%%)", elev_n, 100 * elev_n / lead_n) else "N/A"
  norm_txt <- if (lead_n > 0) sprintf("%d (%.2f%%)", norm_n, 100 * norm_n / lead_n) else "N/A"
  tibble(Count = n, `Lead Elevated (≥5 µg/dL)` = elev_txt, `Lead Normal (<5 µg/dL)` = norm_txt)
}

age_data144 <- long_data %>%
  filter(age_months == 144, !is.na(anemic), !is.na(iron_deficient), !is.na(b12_category))

r_total      <- bind_cols(tibble(`Hb status` = "", `Iron status` = "", `B12 status` = ""), compute_group4(age_data144))

anemic_base  <- age_data144 %>% filter(anemic == TRUE)
r_a1 <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Sufficient",     `B12 status` = "Sufficient"),     compute_group4(anemic_base %>% filter(iron_deficient == FALSE, b12_category == "Sufficient")))
r_a2 <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Sufficient",     `B12 status` = "Not Sufficient"), compute_group4(anemic_base %>% filter(iron_deficient == FALSE, b12_category %in% c("Deficient","Insufficient"))))
r_a3 <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Sufficient"),     compute_group4(anemic_base %>% filter(iron_deficient == TRUE,  b12_category == "Sufficient")))
r_a4 <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Not Sufficient"), compute_group4(anemic_base %>% filter(iron_deficient == TRUE,  b12_category %in% c("Deficient","Insufficient"))))
r_a_total <- bind_cols(tibble(`Hb status` = "Total Anemic", `Iron status` = "", `B12 status` = ""), compute_group4(anemic_base))

nonanemic_base <- age_data144 %>% filter(anemic == FALSE)
r_n1 <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Sufficient",     `B12 status` = "Sufficient"),     compute_group4(nonanemic_base %>% filter(iron_deficient == FALSE, b12_category == "Sufficient")))
r_n2 <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Sufficient",     `B12 status` = "Not Sufficient"), compute_group4(nonanemic_base %>% filter(iron_deficient == FALSE, b12_category %in% c("Deficient","Insufficient"))))
r_n3 <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Sufficient"),     compute_group4(nonanemic_base %>% filter(iron_deficient == TRUE,  b12_category == "Sufficient")))
r_n4 <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Not Sufficient"), compute_group4(nonanemic_base %>% filter(iron_deficient == TRUE,  b12_category %in% c("Deficient","Insufficient"))))
r_n_total <- bind_cols(tibble(`Hb status` = "Total Non-Anemic", `Iron status` = "", `B12 status` = ""), compute_group4(nonanemic_base))

table4 <- bind_rows(r_total, r_a1, r_a2, r_a3, r_a4, r_a_total, r_n1, r_n2, r_n3, r_n4, r_n_total)

table4_caption <- paste0(
  "Table 4: Combined deficiency of iron, vitamin B12 and blood lead exposure among anemic and non-anemic ",
  "(New cut-off) children at 144 months of age. Restricted to children with a concurrent hemoglobin (Hb) ",
  "measurement AND non-missing iron and B12 data simultaneously (hence N differs from Tables 1-3, which ",
  "report iron and B12 completeness separately)."
)

ft4 <- flextable(table4) %>%
  merge_v(j = "Hb status") %>%
  valign(j = "Hb status", valign = "top", part = "body") %>%
  align(j = 1:3, align = "left",   part = "body") %>%
  align(j = 4:6, align = "center", part = "body") %>%
  align(align = "center", part = "header") %>%
  bold(i = c(1, 6, 11), part = "body") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  bold(part = "header") %>%
  add_header_lines(table4_caption) %>%
  autofit()

doc4 <- read_docx() %>%
  body_add_par(table4_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft4)
print(doc4, target = paste0(out_dir, "Table4_Combined_Deficiency_144m.docx"))

# ============================================================
# TABLE 14: Combined iron/B12/lead deficiency at 144 months WITH
# Wilson 95% CI (companion to Table 4). Adds exact CIs to the lead
# columns for every iron x B12 subgroup, most of which have N < 20.
# ============================================================
compute_group4_ci <- function(fd) {
  n      <- nrow(fd)
  lead_n <- sum(!is.na(fd$lead))
  elev_n <- sum(fd$lead >= 5, na.rm = TRUE)
  norm_n <- sum(fd$lead < 5,  na.rm = TRUE)
  elev_txt <- if (lead_n > 0) fmt_pct_ci(elev_n, lead_n) else "N/A"
  norm_txt <- if (lead_n > 0) fmt_pct_ci(norm_n, lead_n) else "N/A"
  tibble(Count = n, `Lead Elevated (≥5 µg/dL)` = elev_txt, `Lead Normal (<5 µg/dL)` = norm_txt)
}

r_total_ci      <- bind_cols(tibble(`Hb status` = "", `Iron status` = "", `B12 status` = ""), compute_group4_ci(age_data144))

r_a1_ci <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Sufficient",     `B12 status` = "Sufficient"),     compute_group4_ci(anemic_base %>% filter(iron_deficient == FALSE, b12_category == "Sufficient")))
r_a2_ci <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Sufficient",     `B12 status` = "Not Sufficient"), compute_group4_ci(anemic_base %>% filter(iron_deficient == FALSE, b12_category %in% c("Deficient","Insufficient"))))
r_a3_ci <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Sufficient"),     compute_group4_ci(anemic_base %>% filter(iron_deficient == TRUE,  b12_category == "Sufficient")))
r_a4_ci <- bind_cols(tibble(`Hb status` = "Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Not Sufficient"), compute_group4_ci(anemic_base %>% filter(iron_deficient == TRUE,  b12_category %in% c("Deficient","Insufficient"))))
r_a_total_ci <- bind_cols(tibble(`Hb status` = "Total Anemic", `Iron status` = "", `B12 status` = ""), compute_group4_ci(anemic_base))

r_n1_ci <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Sufficient",     `B12 status` = "Sufficient"),     compute_group4_ci(nonanemic_base %>% filter(iron_deficient == FALSE, b12_category == "Sufficient")))
r_n2_ci <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Sufficient",     `B12 status` = "Not Sufficient"), compute_group4_ci(nonanemic_base %>% filter(iron_deficient == FALSE, b12_category %in% c("Deficient","Insufficient"))))
r_n3_ci <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Sufficient"),     compute_group4_ci(nonanemic_base %>% filter(iron_deficient == TRUE,  b12_category == "Sufficient")))
r_n4_ci <- bind_cols(tibble(`Hb status` = "Non-Anemic", `Iron status` = "Not Sufficient", `B12 status` = "Not Sufficient"), compute_group4_ci(nonanemic_base %>% filter(iron_deficient == TRUE,  b12_category %in% c("Deficient","Insufficient"))))
r_n_total_ci <- bind_cols(tibble(`Hb status` = "Total Non-Anemic", `Iron status` = "", `B12 status` = ""), compute_group4_ci(nonanemic_base))

table14 <- bind_rows(r_total_ci, r_a1_ci, r_a2_ci, r_a3_ci, r_a4_ci, r_a_total_ci, r_n1_ci, r_n2_ci, r_n3_ci, r_n4_ci, r_n_total_ci)

table14_caption <- paste0(
  "Table 14: Combined deficiency of iron, vitamin B12 and blood lead exposure among anemic and non-anemic ",
  "(New cut-off) children at 144 months of age, with exact Wilson 95% confidence intervals shown in brackets ",
  "for the lead columns. Companion table to Table 4; restricted to children with a concurrent hemoglobin (Hb) ",
  "measurement AND non-missing iron and B12 data simultaneously. This table gives the exact 95% CI for the ",
  "0/17 vs 17/17-type splits referenced in Table 4 (Issue #3), e.g. a 17/17 (100%) cell has a wide interval ",
  "(roughly 80-100%) that should be reported alongside the point estimate."
)

ft14 <- flextable(table14) %>%
  merge_v(j = "Hb status") %>%
  valign(j = "Hb status", valign = "top", part = "body") %>%
  align(j = 1:3, align = "left",   part = "body") %>%
  align(j = 4:6, align = "center", part = "body") %>%
  align(align = "center", part = "header") %>%
  bold(i = c(1, 6, 11), part = "body") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  bold(part = "header") %>%
  fontsize(size = 8, part = "all") %>%
  add_header_lines(table14_caption) %>%
  autofit()

doc14 <- read_docx() %>%
  body_add_par(table14_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft14)
print(doc14, target = paste0(out_dir, "Table14_Combined_Deficiency_144m_WithCI.docx"))

# ============================================================
# Build panel_data (needed for Table 5 AND the Figure below) —
# includes both old (Hb<11) and new (WHO age/sex-specific) classifications
# for ALL Hb timepoints (7,15,24,36,60,84,108,144), all children.
# ============================================================
old_cutoff <- function(hb) ifelse(is.na(hb), NA, hb < 11)

panel_ages <- c(7, 15, 24, 36, 60, 84, 108, 144)
panel_labs <- c("7","15","24","36","60","84","108","144")

panel_data <- data %>%
  mutate(across(starts_with("Hb"), as.numeric)) %>%
  pivot_longer(cols = starts_with("Hb"), names_to = "age_label", values_to = "hb") %>%
  mutate(
    age_months = case_when(
      age_label == "Hb7"   ~ 7,   age_label == "Hb15"  ~ 15,  age_label == "Hb24"  ~ 24,
      age_label == "Hb36"  ~ 36,  age_label == "Hb60"  ~ 60,  age_label == "Hb84"  ~ 84,
      age_label == "Hb108" ~ 108, age_label == "Hb144" ~ 144,
      TRUE ~ NA_real_
    ),
    CAFSEX     = factor(CAFSEX, levels = c(1, 2), labels = c("Boys", "Girls")),
    anemia_old = old_cutoff(hb),
    anemia_new = mapply(hb_cutoff,   hb, age_months, as.character(CAFSEX)),
    severity   = mapply(hb_severity, hb, age_months, as.character(CAFSEX)),
    severity   = factor(severity, levels = c("No Anaemia", "Mild", "Moderate", "Severe"))
  ) %>%
  filter(!is.na(age_months))

# ============================================================
# TABLE 5: Anemia prevalence (Old vs New cut-off) by sex, per age,
#          with sex-comparison p-value (Boys vs Girls) and McNemar
#          p-value (Old vs New cut-off, paired within the same child)
#
# UPDATED (Fisher's exact fallback):
# compute_sex_chi2() now runs chi-square by default but automatically
# falls back to Fisher's exact test whenever any expected cell count
# is < 5 (Boys vs Girls x Anemic/Not is always a 2x2 table, and several
# age groups here have small cell counts, e.g. 60m and 84m where the
# sex-specific proportions are nearly identical, previously giving a
# Yates-corrected chi-square statistic of exactly 0 and a spurious
# p = 1.000). This mirrors the fallback rule already used in
# haz_ttest-adjacent chi-square checks elsewhere in this script and in
# Table 9, so the same small-N rule is applied consistently everywhere
# a 2x2 group comparison is made in this analysis.
# ============================================================
format_p <- function(p) {
  if (is.null(p) || is.na(p)) return("—")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_pct_n <- function(n, total) {
  if (is.na(total) || total == 0) return("—")
  paste0(round(100 * n / total, 1), "% (", n, "/", total, ")")
}

compute_sex_chi2 <- function(df, anemia_col) {
  sub <- df %>% filter(!is.na(.data[[anemia_col]]), !is.na(CAFSEX))
  if (nrow(sub) == 0) return("—")
  tab <- table(
    factor(sub$CAFSEX, levels = c("Boys", "Girls")),
    factor(sub[[anemia_col]], levels = c(FALSE, TRUE))
  )
  res <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
  test_used <- "chi2"
  # Fall back to Fisher's exact test when chi-square's asymptotic
  # assumption is unreliable (any expected cell count < 5). This
  # removes the spurious p = 1.000 values that Yates-corrected
  # chi-square can produce when group proportions are nearly identical.
  if (is.null(res) || any(res$expected < 5)) {
    res <- tryCatch(fisher.test(tab), error = function(e) NULL)
    test_used <- "fisher"
  }
  if (is.null(res)) return("—")
  p_str <- format_p(res$p.value)
  # Append a marker (superscript "f") whenever Fisher's exact test was
  # used instead of chi-square, so the test used is visible directly in
  # the table cell without needing to re-derive it later. See caption
  # for the legend.
  if (test_used == "fisher") paste0(p_str, "\u1da0") else p_str
}

compute_mcnemar <- function(df, col_old, col_new) {
  sub <- df %>% filter(!is.na(.data[[col_old]]), !is.na(.data[[col_new]]))
  if (nrow(sub) == 0) return("—")
  tab <- table(
    factor(sub[[col_old]], levels = c(FALSE, TRUE)),
    factor(sub[[col_new]], levels = c(FALSE, TRUE))
  )
  b <- tab[1,2]; c <- tab[2,1]  # discordant cells
  if ((b + c) < 25) {
    # exact binomial test on discordant pairs, no continuity artifact
    res <- tryCatch(binom.test(min(b,c), b + c, p = 0.5), error = function(e) NULL)
  } else {
    res <- tryCatch(suppressWarnings(mcnemar.test(tab, correct = TRUE)), error = function(e) NULL)
  }
  if (is.null(res)) return("—")
  format_p(res$p.value)
}

build_table5_rows <- function(age_val) {
  age_df   <- panel_data %>% filter(age_months == age_val, !is.na(hb))
  total_n  <- nrow(age_df)
  boys_df  <- age_df %>% filter(CAFSEX == "Boys")
  girls_df <- age_df %>% filter(CAFSEX == "Girls")
  
  old_all_n   <- sum(age_df$anemia_old == TRUE,   na.rm = TRUE); old_all_tot  <- sum(!is.na(age_df$anemia_old))
  old_boy_n   <- sum(boys_df$anemia_old == TRUE,  na.rm = TRUE); old_boy_tot  <- sum(!is.na(boys_df$anemia_old))
  old_girl_n  <- sum(girls_df$anemia_old == TRUE, na.rm = TRUE); old_girl_tot <- sum(!is.na(girls_df$anemia_old))
  row_old <- tibble(
    `Age (Months)` = age_val, Row = "Anemic (Old)",
    `Total N` = as.character(total_n),
    `All children` = fmt_pct_n(old_all_n, old_all_tot),
    Boys  = fmt_pct_n(old_boy_n, old_boy_tot),
    Girls = fmt_pct_n(old_girl_n, old_girl_tot),
    `P val - boys vs girls` = compute_sex_chi2(age_df, "anemia_old")
  )
  
  new_all_n   <- sum(age_df$anemia_new == TRUE,   na.rm = TRUE); new_all_tot  <- sum(!is.na(age_df$anemia_new))
  new_boy_n   <- sum(boys_df$anemia_new == TRUE,  na.rm = TRUE); new_boy_tot  <- sum(!is.na(boys_df$anemia_new))
  new_girl_n  <- sum(girls_df$anemia_new == TRUE, na.rm = TRUE); new_girl_tot <- sum(!is.na(girls_df$anemia_new))
  row_new <- tibble(
    `Age (Months)` = age_val, Row = "Anemic (New)",
    `Total N` = as.character(total_n),
    `All children` = fmt_pct_n(new_all_n, new_all_tot),
    Boys  = fmt_pct_n(new_boy_n, new_boy_tot),
    Girls = fmt_pct_n(new_girl_n, new_girl_tot),
    `P val - boys vs girls` = compute_sex_chi2(age_df, "anemia_new")
  )
  
  row_pval <- tibble(
    `Age (Months)` = age_val, Row = "P val \u2013 old vs new",
    `Total N` = "",
    `All children` = compute_mcnemar(age_df,   "anemia_old", "anemia_new"),
    Boys  = compute_mcnemar(boys_df,  "anemia_old", "anemia_new"),
    Girls = compute_mcnemar(girls_df, "anemia_old", "anemia_new"),
    `P val - boys vs girls` = ""
  )
  
  bind_rows(row_old, row_new, row_pval)
}

table5_ages <- c(7, 15, 24, 36, 60, 84, 108, 144)
table5 <- bind_rows(lapply(table5_ages, build_table5_rows))

table5_caption <- paste0(
  "Table 5: Prevalence of anemia using old (Hb < 11) and new (WHO age/sex-specific) cut-offs, overall and by ",
  "sex, with p-values for the boys vs girls comparison and McNemar p-values (old vs new cut-off) at each age. ",
  "The boys vs girls comparison uses chi-square by default, falling back to Fisher's exact test whenever any ",
  "expected cell count was <5 (this avoids the exact p = 1.000 values that Yates-corrected chi-square can give ",
  "for near-identical group proportions). A superscript \u1da0 marks p-values computed with Fisher's exact test; ",
  "all other boys vs girls p-values are chi-square. Note that p = 1.000 can still occur under Fisher's exact ",
  "test itself when the observed table is the single most probable table given the row/column totals (e.g. the ",
  "60- and 84-month comparisons under the new cut-off) \u2014 this is a genuine exact result, not a computational ",
  "artifact. Note: old and new cut-offs coincide by construction for the 24-59 month band, so values and ",
  "p-values are identical at 24m and 36m."
)

ft5 <- flextable(table5) %>%
  merge_v(j = "Age (Months)") %>%
  valign(j = "Age (Months)", valign = "top", part = "body") %>%
  align(align = "center", part = "all") %>%
  align(j = c("Age (Months)", "Row"), align = "left", part = "body") %>%
  color(i = ~ grepl("^P val", Row), color = "#1F4E96", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption(table5_caption)

doc5 <- read_docx() %>%
  body_add_par(table5_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft5)
print(doc5, target = paste0(out_dir, "Table5_Anemia_OldNew_BySex.docx"))


# ============================================================
# TABLE 10: Anemia prevalence with Wilson 95% CI (old vs new cut-off)
# Standalone — reuses panel_data from the main script above.
# ============================================================
wilson_ci <- function(x, n, z = 1.96) {
  if (is.na(n) || n == 0) return(c(NA, NA))
  p <- x / n
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2))) / denom
  c(lo = max(0, center - margin), hi = min(1, center + margin))
}

fmt_prev_ci <- function(x, n) {
  if (is.na(n) || n == 0) return("—")
  p  <- 100 * x / n
  ci <- 100 * wilson_ci(x, n)
  sprintf("%.1f%% (%d/%d) [%.1f–%.1f]", p, x, n, ci["lo"], ci["hi"])
}

build_table10_rows <- function(age_val) {
  age_df   <- panel_data %>% filter(age_months == age_val, !is.na(hb))
  boys_df  <- age_df %>% filter(CAFSEX == "Boys")
  girls_df <- age_df %>% filter(CAFSEX == "Girls")
  
  row_old <- tibble(
    `Age (Months)` = age_val,
    Cutoff = "Old (Hb < 11)",
    `All children` = fmt_prev_ci(sum(age_df$anemia_old,  na.rm = TRUE), sum(!is.na(age_df$anemia_old))),
    Boys  = fmt_prev_ci(sum(boys_df$anemia_old,  na.rm = TRUE), sum(!is.na(boys_df$anemia_old))),
    Girls = fmt_prev_ci(sum(girls_df$anemia_old, na.rm = TRUE), sum(!is.na(girls_df$anemia_old)))
  )
  
  row_new <- tibble(
    `Age (Months)` = age_val,
    Cutoff = "New (WHO age/sex-specific)",
    `All children` = fmt_prev_ci(sum(age_df$anemia_new,  na.rm = TRUE), sum(!is.na(age_df$anemia_new))),
    Boys  = fmt_prev_ci(sum(boys_df$anemia_new,  na.rm = TRUE), sum(!is.na(boys_df$anemia_new))),
    Girls = fmt_prev_ci(sum(girls_df$anemia_new, na.rm = TRUE), sum(!is.na(girls_df$anemia_new)))
  )
  
  bind_rows(row_old, row_new)
}

table10 <- bind_rows(lapply(panel_ages, build_table10_rows))

ft10 <- flextable(table10) %>%
  merge_v(j = "Age (Months)") %>%
  valign(j = "Age (Months)", valign = "top", part = "body") %>%
  align(align = "center", part = "all") %>%
  align(j = c("Age (Months)", "Cutoff"), align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption("Table 10: Anemia prevalence (old vs. new WHO cut-off), overall and by sex, with Wilson 95% confidence intervals, at each age.")

doc10 <- read_docx() %>%
  body_add_par("Table 10: Anemia prevalence with Wilson 95% CI, old vs. new WHO cut-off, by age and sex.", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft10)
print(doc10, target = paste0(out_dir, "Table10_Anemia_Prevalence_CI.docx"))

# ============================================================
# TABLE 6 / FIGURE: GEE trend model — anemia prevalence over age
# ============================================================
gee_data <- panel_data %>%
  filter(!is.na(pid)) %>%
  arrange(pid, age_months)

fit_gee_linear <- function(outcome) {
  d    <- gee_data %>% filter(!is.na(.data[[outcome]])) %>% arrange(pid, age_months)
  form <- as.formula(paste0(outcome, " ~ age_months"))
  geeglm(form, id = pid, data = d, family = binomial, corstr = "exchangeable")
}

fit_gee_spline <- function(outcome, spline_df = 4) {
  d    <- gee_data %>% filter(!is.na(.data[[outcome]])) %>% arrange(pid, age_months)
  form <- as.formula(paste0(outcome, " ~ ns(age_months, df = ", spline_df, ")"))
  geeglm(form, id = pid, data = d, family = binomial, corstr = "exchangeable")
}

fit_gee_null <- function(outcome) {
  d    <- gee_data %>% filter(!is.na(.data[[outcome]])) %>% arrange(pid, age_months)
  geeglm(as.formula(paste0(outcome, " ~ 1")), id = pid, data = d, family = binomial, corstr = "exchangeable")
}

extract_linear_result <- function(model, label) {
  s     <- summary(model)$coefficients
  est   <- s["age_months", "Estimate"]
  se    <- s["age_months", "Std.err"]
  p     <- s["age_months", "Pr(>|W|)"]
  or_yr <- exp(est * 12)
  lo_yr <- exp((est - 1.96 * se) * 12)
  hi_yr <- exp((est + 1.96 * se) * 12)
  tibble(
    Cutoff = label,
    `OR per year of age` = sprintf("%.3f", or_yr),
    `95% CI` = sprintf("%.3f\u2013%.3f", lo_yr, hi_yr),
    `Linear trend P` = format_p(p)
  )
}

extract_omnibus_p <- function(model_null, model_spline) {
  res <- tryCatch(anova(model_spline, model_null), error = function(e) NULL)
  if (is.null(res)) return("—")
  df_res <- as.data.frame(res)
  pcol   <- names(df_res)[grepl("^p", names(df_res), ignore.case = TRUE)]
  if (length(pcol) == 0) return("—")
  format_p(df_res[[pcol[1]]][1])
}

gee_lin_old    <- fit_gee_linear("anemia_old")
gee_lin_new    <- fit_gee_linear("anemia_new")
gee_spline_old <- fit_gee_spline("anemia_old")
gee_spline_new <- fit_gee_spline("anemia_new")
gee_null_old   <- fit_gee_null("anemia_old")
gee_null_new   <- fit_gee_null("anemia_new")

gee_trend_table <- bind_rows(
  extract_linear_result(gee_lin_old, "Old (Hb < 11)") %>%
    mutate(`Omnibus (nonlinear) age effect P` = extract_omnibus_p(gee_null_old, gee_spline_old)),
  extract_linear_result(gee_lin_new, "New (WHO age/sex-specific)") %>%
    mutate(`Omnibus (nonlinear) age effect P` = extract_omnibus_p(gee_null_new, gee_spline_new))
)

ft6 <- flextable(gee_trend_table) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  align(align = "center", part = "all") %>%
  align(j = 1, align = "left", part = "body") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption("Table 6: GEE logistic regression of anemia on age, accounting for repeated measures within each child (exchangeable working correlation, robust SEs). The linear trend is expressed as an odds ratio per additional year of age; the omnibus P tests the overall (including nonlinear) age effect using a restricted cubic spline (4 df) compared to an intercept-only model.")

doc6 <- read_docx() %>%
  body_add_par("Table 6: GEE logistic regression of anemia on age, accounting for repeated measures within each child.", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft6)
print(doc6, target = paste0(out_dir, "Table6_GEE_Anemia_Trend.docx"))

predict_curve <- function(model, label, age_min = 7, age_max = 144) {
  age_grid <- seq(age_min, age_max, by = 1)
  newdata  <- data.frame(age_months = age_grid)
  
  tt <- delete.response(terms(model))
  mf <- model.frame(tt, newdata, xlev = model$xlevels)
  X  <- model.matrix(tt, mf)
  
  beta <- coef(model)
  V    <- vcov(model)
  
  eta <- as.numeric(X %*% beta)
  se  <- sqrt(pmax(diag(X %*% V %*% t(X)), 0))
  
  tibble(
    age_months = age_grid,
    prevalence = plogis(eta) * 100,
    lo         = plogis(eta - 1.96 * se) * 100,
    hi         = plogis(eta + 1.96 * se) * 100,
    cutoff     = label
  )
}

curve_data <- bind_rows(
  predict_curve(gee_spline_old, "Old (Hb < 11)"),
  predict_curve(gee_spline_new, "New (WHO age/sex-specific)")
)

p_trend <- ggplot(curve_data, aes(x = age_months, y = prevalence, color = cutoff, fill = cutoff)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  scale_color_manual(name = "Cutoff", values = c("New (WHO age/sex-specific)" = "#C2006A", "Old (Hb < 11)" = "#111111")) +
  scale_fill_manual(name = "Cutoff",  values = c("New (WHO age/sex-specific)" = "#C2006A", "Old (Hb < 11)" = "#111111")) +
  scale_y_continuous(limits = c(0, 80), breaks = seq(0, 80, 20), labels = function(x) paste0(x, "%")) +
  scale_x_continuous(breaks = panel_ages) +
  labs(
    title = "Predicted anemia prevalence by age (GEE spline model, with 95% CI)",
    subtitle = "Accounts for repeated measures within child; shaded band = 95% CI",
    x = "Age (months)", y = "Predicted prevalence (%)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    plot.subtitle = element_text(size = 8, color = "grey30"),
    legend.position = "bottom"
  )

ggsave(paste0(out_dir, "Figure_GEE_Anemia_Trend.tiff"),
       plot = p_trend, width = 7, height = 5, dpi = 200, device = "tiff")

# ============================================================
# TABLE 7 / FIGURES: Hb vs. biomarker correlations (Spearman), by age
# ============================================================
biomarker_labels <- c(
  bodyiron = "Body Iron (mg/kg)",
  b12      = "Vitamin B12 (pg/ml)",
  lead     = "Lead (\u00b5g/dL)"
)

age_labels7 <- c(`7` = "7m", `15` = "15m", `24` = "24m", `108` = "108m", `144` = "144m")

biomarker_long <- long_data %>%
  select(pid, age_months, hb, bodyiron, b12, lead) %>%
  pivot_longer(cols = c(bodyiron, b12, lead), names_to = "biomarker", values_to = "value") %>%
  filter(!is.na(hb), !is.na(value)) %>%
  mutate(
    biomarker_label = factor(biomarker_labels[biomarker], levels = biomarker_labels),
    age_label = factor(age_labels7[as.character(age_months)], levels = age_labels7)
  )

spearman_safe <- function(x, y) {
  if (length(x) < 3) return(list(r = NA_real_, p = NA_real_))
  res <- tryCatch(suppressWarnings(cor.test(x, y, method = "spearman")),
                  error = function(e) NULL)
  if (is.null(res)) return(list(r = NA_real_, p = NA_real_))
  list(r = unname(res$estimate), p = res$p.value)
}

cor_results <- biomarker_long %>%
  group_by(age_months, age_label, biomarker, biomarker_label) %>%
  summarise(
    n   = n(),
    res = list(spearman_safe(hb, value)),
    .groups = "drop"
  ) %>%
  mutate(
    rho = map_dbl(res, "r"),
    p   = map_dbl(res, "p")
  ) %>%
  select(-res) %>%
  filter(n >= 3) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  arrange(age_months, biomarker)

get_stars <- function(p) {
  if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""
}

table7 <- cor_results %>%
  mutate(
    Age = age_labels7[as.character(age_months)],
    Biomarker = as.character(biomarker_label),
    N = n,
    `Spearman rho` = sprintf("%.3f", rho),
    `P value` = sapply(p, format_p),
    `FDR-adjusted P` = paste0(sapply(p_adj, format_p), sapply(p_adj, get_stars))
  ) %>%
  select(Age, Biomarker, N, `Spearman rho`, `P value`, `FDR-adjusted P`)

ft7 <- flextable(table7) %>%
  merge_v(j = "Age") %>%
  valign(j = "Age", valign = "top", part = "body") %>%
  align(align = "center", part = "all") %>%
  align(j = c("Age", "Biomarker"), align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption("Table 7: Spearman correlation between hemoglobin and body iron, vitamin B12, and blood lead at each age. Restricted to children with a concurrent Hb measurement. P values are FDR (Benjamini-Hochberg) adjusted across all tests; *** p<0.001, ** p<0.01, * p<0.05.")

doc7 <- read_docx() %>%
  body_add_par("Table 7: Spearman correlation between hemoglobin and body iron, vitamin B12, and blood lead at each age.", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft7)
print(doc7, target = paste0(out_dir, "Table7_Hb_Biomarker_Correlations.docx"))

cor_results_full <- cor_results %>%
  ungroup() %>%
  complete(age_label = levels(biomarker_long$age_label),
           biomarker_label = levels(biomarker_long$biomarker_label)) %>%
  mutate(age_label = factor(age_label, levels = levels(biomarker_long$age_label)),
         biomarker_label = factor(biomarker_label, levels = levels(biomarker_long$biomarker_label)))

p_heatmap <- ggplot(cor_results_full, aes(x = age_label, y = biomarker_label, fill = rho)) +
  geom_tile(color = "grey30", linewidth = 0.8) +
  geom_text(data = filter(cor_results_full, !is.na(rho)),
            aes(label = paste0(sprintf("%.2f", rho), sapply(p_adj, get_stars))),
            size = 3, fontface = "bold") +
  scale_fill_gradient2(name = "Spearman\nrho", low = "#B2182B", mid = "white", high = "#2166AC",
                       midpoint = 0, limits = c(-1, 1), na.value = "grey95") +
  labs(x = "Age", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold")
  )

ggsave(paste0(out_dir, "Figure_Hb_Biomarker_Heatmap.tiff"),
       plot = p_heatmap, width = 6, height = 4, dpi = 200, device = "tiff")

p_scatter <- ggplot(biomarker_long, aes(x = value, y = hb)) +
  geom_point(alpha = 0.35, size = 1, color = "#2A6F97") +
  geom_smooth(method = "lm", se = TRUE, color = "#C2006A", linewidth = 0.6) +
  facet_grid(biomarker_label ~ age_label, scales = "free") +
  labs(title = "Hemoglobin vs. Biomarkers, by Age",
       x = NULL, y = "Hemoglobin (g/dL)") +
  theme_bw(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10),
    strip.text = element_text(face = "bold", size = 7),
    axis.text = element_text(size = 6)
  )

ggsave(paste0(out_dir, "Figure_Hb_Biomarker_Scatter.tiff"),
       plot = p_scatter, width = 10, height = 6, dpi = 200, device = "tiff")

# ============================================================
# TABLE 8: Data completeness — N (%) with data for each variable
# at each age.
#
# UPDATED DENOMINATOR/NUMERATOR CONVENTION (fixes >100% bug):
# Denominator ("Assessed (N)") = N children with a non-missing Hb value
# at that visit, as before.
# Numerator (N for each variable) is now ALSO restricted to children
# who have a concurrent, non-missing Hb value at that same visit —
# i.e., N = count of children with BOTH the biomarker present AND Hb
# present at that age. This guarantees N <= Assessed (N) for every row,
# so no percentage can exceed 100%, and makes the denominators here
# identical in logic to those used to build long_data (Tables 1, 3, 7).
# Children who have (e.g.) Lead/B12/TfR data but a missing Hb at that
# specific visit are excluded from these counts; this affects a small
# number of children (2 at 144 months) whose Hb channel failed or was
# not recorded despite other assays being run.
# ============================================================
hb_col_by_age <- c(`7` = "Hb7", `15` = "Hb15", `24` = "Hb24", `36` = "Hb36",
                   `60` = "Hb60", `84` = "Hb84", `108` = "Hb108", `144` = "Hb144")

assessed_n_by_age <- sapply(hb_col_by_age, function(col) sum(!is.na(data[[col]])))
names(assessed_n_by_age) <- names(hb_col_by_age)

var_map <- tribble(
  ~age_months, ~Age,        ~Variable,                    ~column,
  7,   "7 months",   "Hemoglobin",                 "Hb7",
  7,   "7 months",   "Ferritin",                   "Ferritin7",
  7,   "7 months",   "Transferrin Receptor (TfR)", "tfr7",
  7,   "7 months",   "Body Iron",                  "bodyiron7",
  15,  "15 months",  "Hemoglobin",                 "Hb15",
  15,  "15 months",  "Ferritin",                   "Ferritin15",
  15,  "15 months",  "Transferrin Receptor (TfR)", "tfr15",
  15,  "15 months",  "Body Iron",                  "bodyiron15",
  15,  "15 months",  "Lead",                       "Lead_15",
  24,  "24 months",  "Hemoglobin",                 "Hb24",
  24,  "24 months",  "Ferritin",                   "Ferritin24",
  24,  "24 months",  "Transferrin Receptor (TfR)", "tfr24",
  24,  "24 months",  "Body Iron",                  "bodyiron24",
  24,  "24 months",  "Lead",                       "Lead_24",
  36,  "36 months",  "Hemoglobin",                 "Hb36",
  60,  "60 months",  "Hemoglobin",                 "Hb60",
  84,  "84 months",  "Hemoglobin",                 "Hb84",
  108, "108 months", "Hemoglobin",                 "Hb108",
  108, "108 months", "Ferritin",                   "Ferritin108",
  108, "108 months", "Transferrin Receptor (TfR)", "tfr108",
  108, "108 months", "Body Iron",                  "bodyiron108",
  108, "108 months", "Vitamin B12",                "B12_9y",
  144, "144 months", "Hemoglobin",                 "Hb144",
  144, "144 months", "Ferritin",                   "Ferritin144",
  144, "144 months", "Transferrin Receptor (TfR)", "tfr144",
  144, "144 months", "Body Iron",                  "bodyiron144",
  144, "144 months", "Vitamin B12",                "B12_12y",
  144, "144 months", "Lead",                       "Lead_12y"
) %>%
  filter(column %in% names(data)) %>%
  mutate(hb_column = hb_col_by_age[as.character(age_months)])

completeness_table <- var_map %>%
  rowwise() %>%
  mutate(
    assessed_n = assessed_n_by_age[as.character(age_months)],
    # Restrict the numerator to children with a concurrent (non-missing)
    # Hb measurement at this same age, matching the long_data convention
    # used in Tables 1, 3, and 7. This is the key fix: it guarantees
    # N <= assessed_n, so Pct can never exceed 100%.
    N   = sum(!is.na(data[[column]]) & !is.na(data[[hb_column]])),
    Pct = round(100 * N / assessed_n, 1),
    cell = paste0(N, " (", Pct, "%)")
  ) %>%
  ungroup() %>%
  mutate(Variable = factor(Variable, levels = c("Hemoglobin", "Ferritin", "Transferrin Receptor (TfR)", "Body Iron", "Vitamin B12", "Lead"))) %>%
  select(age_months, Age, assessed_n, Variable, cell) %>%
  arrange(age_months, Variable)

table8 <- completeness_table %>%
  pivot_wider(names_from = Variable, values_from = cell) %>%
  arrange(age_months) %>%
  rename(`Assessed (N)` = assessed_n) %>%
  select(-age_months)

table8_caption <- paste0(
  "Table 8: Data completeness by age. \"Assessed (N)\" is the number of children with a hemoglobin ",
  "measurement at that visit. All other cells show N (%) of children with BOTH that variable AND a ",
  "concurrent hemoglobin measurement, out of the number assessed at that visit \u2014 the same ",
  "concurrent-Hb restriction used to build Tables 1, 3, and 7. A small number of children per visit ",
  "(e.g., 2 at 144 months) have a non-missing biomarker value but a missing Hb at that same visit; ",
  "these children are excluded from the corresponding N here so that percentages cannot exceed 100%. ",
  "Blank cells indicate the variable was not collected at that age."
)

ft8 <- flextable(table8) %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  align(align = "center", part = "all") %>%
  align(j = "Age", align = "left", part = "body") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption(table8_caption)

doc8 <- read_docx() %>%
  body_add_par(table8_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft8)
print(doc8, target = paste0(out_dir, "Table8_Data_Completeness.docx"))

# ============================================================
# TABLE 9: MAR sensitivity check — attrition comparison,
# evaluated separately AT EACH FOLLOW-UP TIMEPOINT (15-144 months).
#
# UPDATED (Fisher's exact fallback):
# safe_chisq_p() now runs chi-square by default but automatically
# falls back to Fisher's exact test whenever any expected cell count
# is < 5. This matters most here because several "N Dropped" groups
# are very small (e.g. 5 children dropped at 15 months, 8 at 24
# months), which is exactly where the chi-square approximation is
# unreliable and previously produced a spurious p = 1.000 (Sex
# p-value at 15 months). The same rule (expected < 5 -> Fisher's) is
# applied consistently in Table 5's compute_sex_chi2().
# ============================================================
followup_ages  <- c(15, 24, 36, 60, 84, 108, 144)
followup_hbcol <- c("Hb15", "Hb24", "Hb36", "Hb60", "Hb84", "Hb108", "Hb144")

attrition_base <- data %>%
  filter(!is.na(Hb7)) %>%
  mutate(
    sex_group = case_when(CAFSEX == 1 ~ "Boys", CAFSEX == 2 ~ "Girls", TRUE ~ NA_character_),
    boys_flag = sex_group == "Boys",
    anemic_7m = mapply(hb_cutoff, Hb7, 7, sex_group)
  )

safe_pct <- function(sub, col) {
  x <- sub[[col]]
  n <- sum(!is.na(x))
  if (n == 0) return("—")
  k <- sum(x == TRUE, na.rm = TRUE)
  sprintf("%.1f%% (%d/%d)", 100 * k / n, k, n)
}

safe_meansd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("—")
  if (length(x) == 1) return(sprintf("%.2f (NA)", x))
  sprintf("%.2f (%.2f)", mean(x), sd(x))
}

safe_ttest_p <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return("—")
  res <- tryCatch(t.test(x, y), error = function(e) NULL)
  if (is.null(res)) return("—")
  format_p(res$p.value)
}

safe_chisq_p <- function(group_flag, outcome_flag) {
  tab <- table(factor(group_flag, levels = c(TRUE, FALSE)), factor(outcome_flag, levels = c(TRUE, FALSE)))
  res <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
  test_used <- "chi2"
  # Fall back to Fisher's exact test when chi-square's asymptotic
  # assumption is unreliable (any expected cell count < 5) — the
  # dropped-group sizes in this table are frequently well below 20,
  # so this fallback is triggered often and is the correct choice
  # here rather than an edge case.
  if (is.null(res) || any(is.na(res$expected)) || any(res$expected < 5)) {
    res <- tryCatch(fisher.test(tab), error = function(e) NULL)
    test_used <- "fisher"
  }
  if (is.null(res)) return("—")
  p_str <- format_p(res$p.value)
  # Append a marker (superscript "f") whenever Fisher's exact test was
  # used instead of chi-square, so the test used is visible directly in
  # the table cell without needing to re-derive it later. See caption
  # for the legend.
  if (test_used == "fisher") paste0(p_str, "\u1da0") else p_str
}

build_timepoint_row <- function(age_val, hb_col) {
  d   <- attrition_base %>% mutate(has_data = !is.na(.data[[hb_col]]))
  ret <- d %>% filter(has_data)
  drp <- d %>% filter(!has_data)
  
  tibble(
    `Follow-up` = paste0(age_val, " months"),
    `N Retained` = nrow(ret),
    `N Dropped`  = nrow(drp),
    `Boys % (Retained)` = safe_pct(ret, "boys_flag"),
    `Boys % (Dropped)`  = safe_pct(drp, "boys_flag"),
    `Sex p-value` = safe_chisq_p(d$has_data, d$boys_flag),
    `Mean Hb7 (Retained)` = safe_meansd(ret$Hb7),
    `Mean Hb7 (Dropped)`  = safe_meansd(drp$Hb7),
    `Hb7 p-value` = safe_ttest_p(ret$Hb7, drp$Hb7),
    `Anemic % 7m (Retained)` = safe_pct(ret, "anemic_7m"),
    `Anemic % 7m (Dropped)`  = safe_pct(drp, "anemic_7m"),
    `Anemia p-value` = safe_chisq_p(d$has_data, d$anemic_7m)
  )
}

table9 <- bind_rows(lapply(seq_along(followup_ages), function(i) build_timepoint_row(followup_ages[i], followup_hbcol[i])))

table9_caption <- paste0(
  "Table 9: MAR sensitivity check - baseline (7-month) sex distribution, mean hemoglobin, and anemia status ",
  "compared between children retained vs. dropped at each individual follow-up timepoint (15-144 months). Sex ",
  "and Anemia p-values use chi-square by default, falling back to Fisher's exact test whenever any expected ",
  "cell count was <5 (this is common here given some dropped-group sizes are as small as 5-8 children, and ",
  "avoids the spurious p = 1.000 that Yates-corrected chi-square can otherwise give for these small tables). ",
  "A superscript \u1da0 marks p-values computed with Fisher's exact test; all other Sex/Anemia p-values are ",
  "chi-square. Note that p = 1.000 can still legitimately occur under Fisher's exact test when the observed ",
  "table is the single most probable table given the row/column totals \u2014 this is an exact result, not a ",
  "computational artifact."
)

ft9 <- flextable(table9) %>%
  align(align = "center", part = "all") %>%
  align(j = "Follow-up", align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 8, part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption(table9_caption)

doc9 <- read_docx() %>%
  body_add_par(table9_caption, style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft9)
print(doc9, target = paste0(out_dir, "Table9_Attrition_Comparison_ByTimepoint.docx"))

# ============================================================
# FIGURE: Anaemia Panel Chart (TIFF) — old vs new WHO cut-offs
# ============================================================
build_panel1 <- function(df, show_x = FALSE) {
  plot_data <- df %>%
    group_by(age_months) %>%
    summarise(
      Total_N      = sum(!is.na(hb)),
      n_anemic_old = sum(anemia_old, na.rm = TRUE),
      n_anemic_new = sum(anemia_new, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      perc_old = round(100 * n_anemic_old / Total_N, 1),
      perc_new = round(100 * n_anemic_new / Total_N, 1)
    ) %>%
    pivot_longer(cols = c(perc_old, perc_new), names_to = "cutoff", values_to = "prevalence") %>%
    mutate(
      cutoff = recode(cutoff, perc_old = "Old (Hb < 11)", perc_new = "New (WHO age-specific)"),
      age_label = factor(age_months, levels = panel_ages, labels = panel_labs)
    )
  
  ggplot(plot_data, aes(x = age_label, y = prevalence, color = cutoff, group = cutoff, shape = cutoff)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = paste0(prevalence, "%")), size = 2.3, fontface = "bold", show.legend = FALSE,
                    box.padding = 0.3, point.padding = 0.2, min.segment.length = 0.1, segment.size = 0.3,
                    segment.color = "grey60", max.overlaps = Inf, direction = "y", nudge_y = 3) +
    scale_color_manual(name = "Cutoff", values = c("New (WHO age-specific)" = "#C2006A", "Old (Hb < 11)" = "#111111")) +
    scale_shape_manual(name = "Cutoff", values = c("New (WHO age-specific)" = 16, "Old (Hb < 11)" = 16)) +
    scale_y_continuous(limits = c(0, 80), breaks = seq(0, 60, 20), labels = function(x) paste0(x, "%")) +
    labs(title = "2.a - Prevalence of Anaemia using old and new WHO cut-offs for All Children",
         x = if (show_x) "Age (months)" else NULL, y = "Prevalence (%)") +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 9),
      axis.title.y = element_text(face = "bold", size = 7),
      axis.title.x = element_text(face = "bold", size = 7),
      axis.text = element_text(face = "bold", size = 7),
      axis.text.x = if (show_x) element_text(face = "bold", size = 7) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 7),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.35, "cm"),
      legend.margin = margin(t = 0, b = 0),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 3, r = 6, b = 0, l = 6)
    )
}

build_panel_sex <- function(df, cutoff_type = c("old", "new"), panel_title, show_x = FALSE) {
  cutoff_type <- match.arg(cutoff_type)
  anemia_col  <- if (cutoff_type == "old") "anemia_old" else "anemia_new"
  
  plot_data <- df %>%
    filter(!is.na(CAFSEX)) %>%
    group_by(age_months, CAFSEX) %>%
    summarise(Total_N = sum(!is.na(hb)), n_anemic = sum(.data[[anemia_col]], na.rm = TRUE), .groups = "drop") %>%
    mutate(
      prevalence = round(100 * n_anemic / Total_N, 1),
      age_label  = factor(age_months, levels = panel_ages, labels = panel_labs)
    )
  
  ggplot(plot_data, aes(x = age_label, y = prevalence, color = CAFSEX, group = CAFSEX, shape = CAFSEX)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    geom_text_repel(aes(label = paste0(prevalence, "%")), size = 2.3, fontface = "bold", show.legend = FALSE,
                    box.padding = 0.3, point.padding = 0.2, min.segment.length = 0.1, segment.size = 0.3,
                    segment.color = "grey60", max.overlaps = Inf, direction = "y", nudge_y = 3) +
    scale_color_manual(name = "Sex", values = c("Boys" = "#111111", "Girls" = "#C2006A")) +
    scale_shape_manual(name = "Sex", values = c("Boys" = 16, "Girls" = 17)) +
    scale_y_continuous(limits = c(0, 80), breaks = seq(0, 60, 20), labels = function(x) paste0(x, "%")) +
    labs(title = panel_title, x = if (show_x) "Age (months)" else NULL, y = "Prevalence (%)") +
    theme_bw(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0, size = 9),
      axis.title.y = element_text(face = "bold", size = 7),
      axis.title.x = element_text(face = "bold", size = 7),
      axis.text = element_text(face = "bold", size = 7),
      axis.text.x = if (show_x) element_text(face = "bold", size = 7) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 7),
      legend.text = element_text(size = 7),
      legend.key.size = unit(0.35, "cm"),
      legend.margin = margin(t = 0, b = 0),
      panel.grid.minor = element_blank(),
      plot.margin = margin(t = 3, r = 6, b = 0, l = 6)
    )
}

p1 <- build_panel1(panel_data, show_x = FALSE)
p2 <- build_panel_sex(panel_data, cutoff_type = "old",
                      panel_title = "2.b - Prevalence of Anaemia using old cut-off (Hb < 11) by sex", show_x = FALSE)
p3 <- build_panel_sex(panel_data, cutoff_type = "new",
                      panel_title = "2.c - Prevalence of Anaemia using New WHO cut-off (age & sex-specific) by sex", show_x = TRUE)

final_panel <- (p1 / p2 / p3) + plot_layout(heights = c(1, 1, 1))

ggsave(paste0(out_dir, "Figure_Anemia_Panel_Chart.tiff"),
       plot = final_panel, width = 6, height = 7, dpi = 200, device = "tiff")

message("All outputs saved to: ", out_dir)

# ============================================================
# TABLE 11: Nutritional status (HAZ) by anemia status
#
# NOTE: This table reports ONLY the independent-samples t-test
# comparing mean HAZ between Anemic and Non-Anemic children at each
# age (see haz_ttest_p() below and the "P value" column in table11).
# There is no categorical (stunted vs not) comparison here, so no
# chi-square/Fisher's exact test is needed or used in this table —
# any earlier chi-square helper for HAZ has been removed as unused.
# ============================================================
haz_cols <- c(`7` = "HAZ_7", `15` = "HAZ_15", `24` = "HAZ_24",
              `108` = "HAZ_108", `144` = "HAZ_144")

haz_long <- data %>%
  select(pid, all_of(unname(haz_cols))) %>%
  pivot_longer(cols = -pid, names_to = "haz_label", values_to = "haz") %>%
  mutate(
    age_months = as.numeric(names(haz_cols)[match(haz_label, haz_cols)]),
    haz = as.numeric(as.character(haz))
  ) %>%
  filter(!is.na(haz)) %>%
  select(pid, age_months, haz)

haz_anemia <- long_data %>%
  select(pid, age_months, hb_status) %>%
  filter(!is.na(hb_status)) %>%
  inner_join(haz_long, by = c("pid", "age_months")) %>%
  mutate(
    stunted        = haz < -2,
    severe_stunted = haz < -3
  )

safe_meansd_haz <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return("—")
  sprintf("%.2f (%.2f)", mean(x), sd(x))
}

haz_ttest_p <- function(df) {
  a <- df$haz[df$hb_status == "Anemic"]
  b <- df$haz[df$hb_status == "Non-Anemic"]
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (length(a) < 2 || length(b) < 2) return("—")
  res <- tryCatch(t.test(a, b), error = function(e) NULL)
  if (is.null(res)) return("—")
  format_p(res$p.value)
}

build_table11_rows <- function(age_val) {
  age_df <- haz_anemia %>% filter(age_months == age_val)
  
  bind_rows(
    tibble(
      `Age (Months)` = age_val, `Hb Status` = "Anemic",
      N = sum(age_df$hb_status == "Anemic"),
      `Mean HAZ (SD)` = safe_meansd_haz(age_df$haz[age_df$hb_status == "Anemic"]),
      `P value` = ""
    ),
    tibble(
      `Age (Months)` = age_val, `Hb Status` = "Non-Anemic",
      N = sum(age_df$hb_status == "Non-Anemic"),
      `Mean HAZ (SD)` = safe_meansd_haz(age_df$haz[age_df$hb_status == "Non-Anemic"]),
      `P value` = haz_ttest_p(age_df)
    )
  )
}

table11 <- bind_rows(lapply(c(7, 15, 24, 108, 144), build_table11_rows)) %>%
  mutate(N = ifelse(is.na(N), "", as.character(N)))

ft11 <- flextable(table11) %>%
  merge_v(j = c("Age (Months)", "P value")) %>%
  valign(j = c("Age (Months)", "P value"), valign = "center", part = "body") %>%
  align(align = "center", part = "all") %>%
  align(j = c("Age (Months)", "Hb Status"), align = "left", part = "body") %>%
  bold(part = "header") %>%
  fontsize(size = 9, part = "all") %>%
  border_outer(part = "all", border = fp_border(color = "black", width = 1)) %>%
  border_inner(part = "all", border = fp_border(color = "black", width = 0.5)) %>%
  autofit() %>%
  set_caption("Table 11: Nutritional status (height-for-age z-score, HAZ) by anemia status (New WHO cut-off) at 7, 15, 24, 108, and 144 months. P value from independent-samples t-test comparing mean HAZ between Anemic and Non-Anemic children.")

doc11 <- read_docx() %>%
  body_add_par("Table 11: Nutritional status (HAZ) by anemia status.", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft11)
print(doc11, target = paste0(out_dir, "Table11_HAZ_by_AnemiaStatus.docx"))














# ============================================================
# DIAGNOSTIC: verify the p = 1.000 cells in Table 5 (boys vs girls)
# and Table 9 (retained vs dropped, sex).
#
# Run this AFTER your main script has run (so panel_data,
# attrition_base, etc. already exist in your environment).
# It prints, for each flagged cell:
#   - the raw 2x2 table
#   - chi-square expected counts (to check the <5 rule)
#   - chi-square p-value (with and without continuity correction)
#   - Fisher's exact p-value
#   - which one your pipeline actually used and why
# ============================================================

inspect_sex_comparison <- function(age_val, anemia_col, label) {
  age_df <- panel_data %>% filter(age_months == age_val, !is.na(hb),
                                  !is.na(.data[[anemia_col]]), !is.na(CAFSEX))
  tab <- table(
    factor(age_df$CAFSEX, levels = c("Boys", "Girls")),
    factor(age_df[[anemia_col]], levels = c(FALSE, TRUE))
  )
  cat("\n==============================\n")
  cat(label, "\n")
  cat("2x2 table (rows = Boys/Girls, cols = Not Anemic/Anemic):\n")
  print(tab)
  
  chi_nc <- suppressWarnings(chisq.test(tab, correct = FALSE))
  chi_c  <- suppressWarnings(chisq.test(tab, correct = TRUE))
  fis    <- fisher.test(tab)
  
  cat("\nExpected counts (chi-square):\n")
  print(round(chi_nc$expected, 2))
  cat("\nAny expected cell < 5? ", any(chi_nc$expected < 5), "\n")
  
  cat(sprintf("\nChi-square (no continuity correction): stat = %.4f, p = %.4f\n",
              chi_nc$statistic, chi_nc$p.value))
  cat(sprintf("Chi-square (Yates continuity correction): stat = %.4f, p = %.4f\n",
              chi_c$statistic, chi_c$p.value))
  cat(sprintf("Fisher's exact test: p = %.4f\n", fis$p.value))
  cat(sprintf("Odds ratio (Fisher): %.3f  [95%% CI: %.3f - %.3f]\n",
              fis$estimate, fis$conf.int[1], fis$conf.int[2]))
}

inspect_retained_dropped <- function(age_val, hb_col, label) {
  d   <- attrition_base %>% mutate(has_data = !is.na(.data[[hb_col]]))
  tab <- table(
    factor(d$has_data,  levels = c(TRUE, FALSE)),
    factor(d$boys_flag, levels = c(TRUE, FALSE))
  )
  cat("\n==============================\n")
  cat(label, "\n")
  cat("2x2 table (rows = Retained/Dropped, cols = Boys/Girls):\n")
  print(tab)
  
  chi_nc <- suppressWarnings(chisq.test(tab, correct = FALSE))
  chi_c  <- suppressWarnings(chisq.test(tab, correct = TRUE))
  fis    <- fisher.test(tab)
  
  cat("\nExpected counts (chi-square):\n")
  print(round(chi_nc$expected, 2))
  cat("\nAny expected cell < 5? ", any(chi_nc$expected < 5), "\n")
  
  cat(sprintf("\nChi-square (no continuity correction): stat = %.4f, p = %.4f\n",
              chi_nc$statistic, chi_nc$p.value))
  cat(sprintf("Chi-square (Yates continuity correction): stat = %.4f, p = %.4f\n",
              chi_c$statistic, chi_c$p.value))
  cat(sprintf("Fisher's exact test: p = %.4f\n", fis$p.value))
}

# ---- The specific cells flagged as p = 1.000 in your output ----

inspect_sex_comparison(60,  "anemia_new", "Table 5 | 60 months | Anemic (New) | Boys vs Girls")
inspect_sex_comparison(84,  "anemia_new", "Table 5 | 84 months | Anemic (New) | Boys vs Girls")
inspect_sex_comparison(108, "anemia_old", "Table 5 | 108 months | Anemic (Old) | Boys vs Girls")

inspect_retained_dropped(15, "Hb15", "Table 9 | 15 months | Sex: Retained vs Dropped")

cat("\n\n==============================\n")
cat("INTERPRETATION\n")
cat("==============================\n")
cat("If 'Any expected cell < 5?' is TRUE, your pipeline correctly used\n")
cat("Fisher's exact test (marked with the superscript in your tables).\n")
cat("If Fisher's p and the uncorrected chi-square p are BOTH close to 1,\n")
cat("that confirms the groups genuinely have near-identical proportions\n")
cat("(or, for small dropped groups, the data are consistent with no\n")
cat("association) -- i.e. p = 1.000 is a real result, not a bug.\n")