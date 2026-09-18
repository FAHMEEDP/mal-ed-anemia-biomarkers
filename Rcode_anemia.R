library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(readxl)
library(geepack)
library(splines)

data <- read_excel("maled_data.xlsx") %>%
  as.data.frame()

data <- data %>%
  rename(
    pid = si_no, CAFSEX = sex,
    Hb7 = Hb_7, Hb15 = Hb_15, Hb24 = Hb_24, Hb36 = Hb_36,
    Hb60 = Hb_60, Hb84 = Hb_84, Hb108 = Hb_108, Hb144 = Hb_144,
    Ferritin7 = Ferritin_7, Ferritin15 = Ferritin_15, Ferritin24 = Ferritin_24,
    Ferritin108 = Ferritin_108, Ferritin144 = Ferritin_144,
    tfr7 = Tfr_7, tfr15 = Tfr_15, tfr24 = Tfr_24, tfr108 = Tfr_108, tfr144 = Tfr_144,
    bodyiron7 = Bodyiron_7, bodyiron15 = Bodyiron_15, bodyiron24 = Bodyiron_24,
    bodyiron108 = Bodyiron_108, bodyiron144 = Bodyiron_144,
    B12_9y = B12_108, B12_12y = B12_144, Lead_12y = Lead_144
  )

num_cols <- c("Hb7","Hb15","Hb24","Hb36","Hb60","Hb84","Hb108","Hb144",
              "Ferritin7","Ferritin15","Ferritin24","Ferritin108","Ferritin144",
              "tfr7","tfr15","tfr24","tfr108","tfr144",
              "bodyiron7","bodyiron15","bodyiron24","bodyiron108","bodyiron144",
              "B12_9y","B12_12y","Lead_15","Lead_24","Lead_12y")
data[num_cols] <- lapply(data[num_cols], function(x) as.numeric(as.character(x)))






# Classification functions - WHO 2024 age/sex-specific Hb cut-offs

hb_cutoff <- function(hb, age_months, sex = NA) {
  if (is.na(hb)) return(NA)
  if (age_months <= 23)  return(hb < 10.5)
  if (age_months <= 59)  return(hb < 11.0)
  if (age_months <= 131) return(hb < 11.5)
  if (age_months <= 167) return(hb < 12.0)
  if (!is.na(sex) && sex == "Boys") return(hb < 13.0) else return(hb < 12.0)
}

hb_severity <- function(hb, age_months, sex = NA) {
  if (is.na(hb)) return(NA_character_)
  if (age_months <= 23)       cuts <- c(10.5, 9.5, 7.0)
  else if (age_months <= 59)  cuts <- c(11.0, 10.0, 7.0)
  else if (age_months <= 131) cuts <- c(11.5, 11.0, 8.0)
  else if (age_months <= 167) cuts <- c(12.0, 11.0, 8.0)
  else cuts <- c(if (!is.na(sex) && sex == "Boys") 13.0 else 12.0, 11.0, 8.0)
  if (hb >= cuts[1]) "No Anaemia" else if (hb >= cuts[2]) "Mild" else if (hb >= cuts[3]) "Moderate" else "Severe"
}

iron_cutoff <- function(iron) ifelse(is.na(iron), NA, iron < 0)
b12_status  <- function(b12) {
  if (is.na(b12)) NA else if (b12 < 200) "Deficient" else if (b12 < 300) "Insufficient" else "Sufficient"
}




## Formatting / statistical-test helpers

wilson_ci <- function(x, n, z = 1.96) {
  if (is.na(n) || n == 0) return(c(lo = NA_real_, hi = NA_real_))
  p <- x / n
  denom  <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  margin <- z * sqrt((p * (1 - p) / n) + (z^2 / (4 * n^2))) / denom
  c(lo = max(0, centre - margin), hi = min(1, centre + margin))
}

fmt_pct <- Vectorize(function(n, total) {
  if (is.na(total) || total == 0) "-" else sprintf("%.1f%% (%d/%d)", 100 * n / total, n, total)
})

fmt_pct_ci <- Vectorize(function(n, total) {
  if (is.na(total) || total == 0) return("-")
  ci <- 100 * wilson_ci(n, total)
  sprintf("%.1f%% (%d/%d) [%.1f-%.1f]", 100 * n / total, n, total, ci["lo"], ci["hi"])
})

format_p <- function(p) if (is.null(p) || is.na(p)) "-" else if (p < 0.001) "<0.001" else sprintf("%.3f", p)


# chi-square with Fisher's-exact fallback whenever any expected cell < 5;
# Fisher's-exact results are marked with a superscript "f" so the test

test_2x2 <- function(tab) {
  chi <- tryCatch(suppressWarnings(chisq.test(tab)), error = function(e) NULL)
  if (is.null(chi) || any(is.na(chi$expected)) || any(chi$expected < 5)) {
    p <- fisher.test(tab)$p.value
    paste0(format_p(p), "\u1da0")
  } else {
    format_p(chi$p.value)
  }
}

sex_test <- function(df, outcome_col) {
  sub <- df %>% filter(!is.na(.data[[outcome_col]]), !is.na(CAFSEX))
  if (nrow(sub) == 0) return("-")
  tab <- table(factor(sub$CAFSEX, levels = c("Boys","Girls")),
               factor(sub[[outcome_col]], levels = c(FALSE, TRUE)))
  test_2x2(tab)
}

attrition_test <- function(has_data, outcome) {
  tab <- table(factor(has_data, levels = c(TRUE, FALSE)), factor(outcome, levels = c(TRUE, FALSE)))
  test_2x2(tab)
}

safe_meansd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2) return("-")
  sprintf("%.2f (%.2f)", mean(x), sd(x))
}

safe_ttest_p <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (length(x) < 2 || length(y) < 2) return("-")
  format_p(t.test(x, y)$p.value)
}



# Comparing prevalence between old and new Hb cut-off (7 & 15 months)

old_band_cutoff <- 11.0  # previous (2011) WHO cut-off for the 6-23 month band
for (age in c(7, 15)) {
  hb_vals <- data[[paste0("Hb", age)]]
  hb_vals <- hb_vals[!is.na(hb_vals)]
  n_total <- length(hb_vals)
  n_anemic_old <- sum(hb_vals < old_band_cutoff)
  cat(sprintf("%d months - previous WHO cut-off (Hb < %.1f g/dL): %.1f%% (%d/%d)\n",
              age, old_band_cutoff, 100 * n_anemic_old / n_total, n_anemic_old, n_total))
}






# Hb-based dataset - all ages, all children

panel_data <- data %>%
  pivot_longer(starts_with("Hb"), names_to = "age_label", values_to = "hb") %>%
  mutate(
    age_months = as.numeric(sub("Hb", "", age_label)),
    CAFSEX = factor(CAFSEX, levels = c(1,2), labels = c("Boys","Girls")),
    anemic = mapply(hb_cutoff, hb, age_months, as.character(CAFSEX)),
    severity = factor(mapply(hb_severity, hb, age_months, as.character(CAFSEX)),
                      levels = c("No Anaemia","Mild","Moderate","Severe"))
  ) %>%
  filter(!is.na(hb))

panel_ages <- c(7, 15, 24, 36, 60, 84, 108, 144)



# Table 1: anemia prevalence + severity by sex, with Wilson 95% CI

build_table1_row <- function(age_val) {
  d     <- panel_data %>% filter(age_months == age_val)
  boys  <- d %>% filter(CAFSEX == "Boys")
  girls <- d %>% filter(CAFSEX == "Girls")
  N <- nrow(d)
  tibble(
    `Age (Months)` = age_val,
    `Total N` = N,
    `All children` = fmt_pct_ci(sum(d$anemic, na.rm = TRUE), N),
    Mild     = fmt_pct_ci(sum(d$severity == "Mild",     na.rm = TRUE), N),
    Moderate = fmt_pct_ci(sum(d$severity == "Moderate", na.rm = TRUE), N),
    Severe   = fmt_pct_ci(sum(d$severity == "Severe",   na.rm = TRUE), N),
    Boys  = fmt_pct(sum(boys$anemic,  na.rm = TRUE), nrow(boys)),
    Girls = fmt_pct(sum(girls$anemic, na.rm = TRUE), nrow(girls)),
    `P val - boys vs girls` = sex_test(d, "anemic")
  )
}

table1 <- bind_rows(lapply(panel_ages, build_table1_row))
print(table1)





# GEE trend of anemia prevalence over age

gee_data <- panel_data %>% filter(!is.na(pid)) %>% arrange(pid, age_months)

fit_gee <- function(spline = FALSE) {
  rhs <- if (spline) "ns(age_months, df = 4)" else "age_months"
  geeglm(as.formula(paste("anemic ~", rhs)), id = pid, data = gee_data,
         family = binomial, corstr = "exchangeable")
}

gee_lin    <- fit_gee()
gee_spline <- fit_gee(spline = TRUE)
gee_null   <- geeglm(anemic ~ 1, id = pid, data = gee_data, family = binomial, corstr = "exchangeable")

coefs <- summary(gee_lin)$coefficients
est <- coefs["age_months", "Estimate"]; se <- coefs["age_months", "Std.err"]
or  <- exp(est * 12); ci <- exp((est + c(-1, 1) * 1.96 * se) * 12)

anova_df <- as.data.frame(anova(gee_spline, gee_null))
omni_p   <- format_p(anova_df[1, grep("^p", names(anova_df), ignore.case = TRUE)])

gee_table <- tibble(
  Cutoff = "WHO age/sex-specific (2024)",
  `OR per year` = sprintf("%.3f", or),
  `95% CI` = sprintf("%.3f-%.3f", ci[1], ci[2]),
  `Linear trend P` = format_p(coefs["age_months", "Pr(>|W|)"]),
  `Omnibus (nonlinear) P` = omni_p
)
print(gee_table)




# Figure 2: anemia prevalence by age, overall and by sex (WHO cut-off only)

fig2a <- panel_data %>%
  group_by(age_months) %>%
  summarise(N = sum(!is.na(hb)), n_anemic = sum(anemic, na.rm = TRUE), .groups = "drop") %>%
  mutate(Prevalence = round(100 * n_anemic / N, 1), age_label = factor(age_months, levels = panel_ages)) %>%
  ggplot(aes(age_label, Prevalence, group = 1)) +
  geom_line(linewidth = 0.9, colour = "#C2006A") + geom_point(size = 2, colour = "#C2006A") +
  geom_text_repel(aes(label = paste0(Prevalence, "%")), size = 2.3, show.legend = FALSE) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(title = "2.a - Prevalence of Anaemia (WHO age/sex-specific cut-off) for All Children",
       x = NULL, y = "Prevalence (%)") +
  theme_bw(base_size = 8)

fig2b <- panel_data %>%
  filter(!is.na(CAFSEX)) %>%
  group_by(age_months, CAFSEX) %>%
  summarise(N = sum(!is.na(hb)), n = sum(anemic, na.rm = TRUE), .groups = "drop") %>%
  mutate(Prevalence = round(100 * n / N, 1), age_label = factor(age_months, levels = panel_ages)) %>%
  ggplot(aes(age_label, Prevalence, colour = CAFSEX, group = CAFSEX)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  geom_text_repel(aes(label = paste0(Prevalence, "%")), size = 2.3, show.legend = FALSE) +
  scale_colour_manual(values = c(Boys = "#111111", Girls = "#C2006A")) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(title = "2.b - Prevalence of Anaemia (WHO age/sex-specific cut-off) by sex",
       x = "Age (months)", y = "Prevalence (%)") +
  theme_bw(base_size = 8) + theme(legend.position = "bottom")

fig2 <- fig2a / fig2b
print(fig2)






# Iron, B12, and lead prevalence are restricted to children with a hemoglobin value at the SAME visit, so denominators match Table 1/Fig 2.

long_data <- data %>%
  select(pid, CAFSEX, Hb7, Hb15, Hb24, Hb108, Hb144,
         bodyiron7, bodyiron15, bodyiron24, bodyiron108, bodyiron144,
         B12_9y, B12_12y, Lead_15, Lead_24, Lead_12y) %>%
  mutate(sex_group = case_when(CAFSEX == 1 ~ "Boys", CAFSEX == 2 ~ "Girls", TRUE ~ NA_character_)) %>%
  pivot_longer(c(Hb7, Hb15, Hb24, Hb108, Hb144), names_to = "age_label", values_to = "hb") %>%
  mutate(age_months = as.numeric(sub("Hb", "", age_label))) %>%
  filter(!is.na(hb)) %>%
  mutate(
    bodyiron = case_when(age_months == 7 ~ bodyiron7, age_months == 15 ~ bodyiron15,
                         age_months == 24 ~ bodyiron24, age_months == 108 ~ bodyiron108,
                         age_months == 144 ~ bodyiron144),
    b12  = case_when(age_months == 108 ~ B12_9y, age_months == 144 ~ B12_12y, TRUE ~ NA_real_),
    lead = case_when(age_months == 15 ~ Lead_15, age_months == 24 ~ Lead_24,
                     age_months == 144 ~ Lead_12y, TRUE ~ NA_real_),
    anemic = mapply(hb_cutoff, hb, age_months, sex_group),
    hb_status = ifelse(anemic, "Anemic", "Non-Anemic"),
    iron_deficient = iron_cutoff(bodyiron),
    b12_category = sapply(b12, b12_status)
  ) %>%
  select(pid, sex_group, age_months, hb, anemic, hb_status, bodyiron, iron_deficient, b12, b12_category, lead)

biomarker_ages <- c(7, 15, 24, 108, 144)






# Table 2: iron / B12 / lead prevalence by anemia status, with Wilson 95% CI

biomarker_row_ci <- function(sub) {
  iron_n <- sum(!is.na(sub$iron_deficient))
  b12_n  <- sum(!is.na(sub$b12_category))
  lead_n <- sum(!is.na(sub$lead))
  tibble(
    `Iron Deficient`   = if (iron_n) fmt_pct_ci(sum(sub$iron_deficient, na.rm = TRUE), iron_n) else "-",
    `Iron Sufficient`  = if (iron_n) fmt_pct_ci(sum(!sub$iron_deficient, na.rm = TRUE), iron_n) else "-",
    `B12 Deficient`    = if (b12_n) fmt_pct_ci(sum(sub$b12_category == "Deficient",    na.rm = TRUE), b12_n) else "-",
    `B12 Insufficient` = if (b12_n) fmt_pct_ci(sum(sub$b12_category == "Insufficient", na.rm = TRUE), b12_n) else "-",
    `B12 Sufficient`   = if (b12_n) fmt_pct_ci(sum(sub$b12_category == "Sufficient",   na.rm = TRUE), b12_n) else "-",
    `Lead Normal`      = if (lead_n) fmt_pct_ci(sum(sub$lead < 5,  na.rm = TRUE), lead_n) else "-",
    `Lead Elevated`    = if (lead_n) fmt_pct_ci(sum(sub$lead >= 5, na.rm = TRUE), lead_n) else "-"
  )
}

table2 <- lapply(biomarker_ages, function(age) {
  bind_rows(
    tibble(Age = age, `Hb Status` = "Anemic")    %>% bind_cols(biomarker_row_ci(long_data %>% filter(age_months == age, hb_status == "Anemic"))),
    tibble(Age = age, `Hb Status` = "Non-Anemic") %>% bind_cols(biomarker_row_ci(long_data %>% filter(age_months == age, hb_status == "Non-Anemic")))
  )
}) %>% bind_rows()
print(table2)





# Supplementary Table 2: iron / B12 / lead prevalence, all children, by age & sex

biomarker_row_pct <- function(sub) {
  iron_n <- sum(!is.na(sub$iron_deficient)); b12_n <- sum(!is.na(sub$b12_category)); lead_n <- sum(!is.na(sub$lead))
  tibble(
    `Iron Deficient`   = if (iron_n) fmt_pct(sum(sub$iron_deficient, na.rm = TRUE), iron_n) else "-",
    `Iron Sufficient`  = if (iron_n) fmt_pct(sum(!sub$iron_deficient, na.rm = TRUE), iron_n) else "-",
    `B12 Deficient`    = if (b12_n) fmt_pct(sum(sub$b12_category == "Deficient",    na.rm = TRUE), b12_n) else "-",
    `B12 Insufficient` = if (b12_n) fmt_pct(sum(sub$b12_category == "Insufficient", na.rm = TRUE), b12_n) else "-",
    `B12 Sufficient`   = if (b12_n) fmt_pct(sum(sub$b12_category == "Sufficient",   na.rm = TRUE), b12_n) else "-",
    `Lead Normal`      = if (lead_n) fmt_pct(sum(sub$lead < 5,  na.rm = TRUE), lead_n) else "-",
    `Lead Elevated`    = if (lead_n) fmt_pct(sum(sub$lead >= 5, na.rm = TRUE), lead_n) else "-"
  )
}

supp_table2 <- lapply(biomarker_ages, function(age) {
  age_df <- long_data %>% filter(age_months == age)
  bind_rows(
    tibble(Age = age, Sex = "All Children") %>% bind_cols(biomarker_row_pct(age_df)),
    tibble(Age = age, Sex = "Boys")  %>% bind_cols(biomarker_row_pct(age_df %>% filter(sex_group == "Boys"))),
    tibble(Age = age, Sex = "Girls") %>% bind_cols(biomarker_row_pct(age_df %>% filter(sex_group == "Girls")))
  )
}) %>% bind_rows()
print(supp_table2)







# Supplementary Table 3: combined iron/B12/lead status at 144 months, with Wilson 95% CI

d144 <- long_data %>% filter(age_months == 144, !is.na(anemic), !is.na(iron_deficient), !is.na(b12_category))

lead_row_ci <- function(sub) {
  n <- sum(!is.na(sub$lead))
  tibble(N = nrow(sub),
         `Lead Elevated` = if (n) fmt_pct_ci(sum(sub$lead >= 5, na.rm = TRUE), n) else "-",
         `Lead Normal`   = if (n) fmt_pct_ci(sum(sub$lead < 5,  na.rm = TRUE), n) else "-")
}

subgroup <- function(anemic_flag, iron_def, b12_lab) {
  d144 %>% filter(anemic == anemic_flag, iron_deficient == iron_def,
                  if (b12_lab == "Sufficient") b12_category == "Sufficient" else b12_category %in% c("Deficient","Insufficient"))
}

supp_table3 <- bind_rows(
  tibble(`Hb status` = "Total", `Iron status` = "", `B12 status` = "") %>% bind_cols(lead_row_ci(d144)),
  tibble(`Hb status` = "Anemic", `Iron status` = "Sufficient", `B12 status` = "Sufficient")         %>% bind_cols(lead_row_ci(subgroup(TRUE,  FALSE, "Sufficient"))),
  tibble(`Hb status` = "Anemic", `Iron status` = "Sufficient", `B12 status` = "Not sufficient")     %>% bind_cols(lead_row_ci(subgroup(TRUE,  FALSE, "Not sufficient"))),
  tibble(`Hb status` = "Anemic", `Iron status` = "Not sufficient", `B12 status` = "Sufficient")     %>% bind_cols(lead_row_ci(subgroup(TRUE,  TRUE,  "Sufficient"))),
  tibble(`Hb status` = "Anemic", `Iron status` = "Not sufficient", `B12 status` = "Not sufficient") %>% bind_cols(lead_row_ci(subgroup(TRUE,  TRUE,  "Not sufficient"))),
  tibble(`Hb status` = "Total Anemic", `Iron status` = "", `B12 status` = "") %>% bind_cols(lead_row_ci(d144 %>% filter(anemic == TRUE))),
  tibble(`Hb status` = "Non-Anemic", `Iron status` = "Sufficient", `B12 status` = "Sufficient")         %>% bind_cols(lead_row_ci(subgroup(FALSE, FALSE, "Sufficient"))),
  tibble(`Hb status` = "Non-Anemic", `Iron status` = "Sufficient", `B12 status` = "Not sufficient")     %>% bind_cols(lead_row_ci(subgroup(FALSE, FALSE, "Not sufficient"))),
  tibble(`Hb status` = "Non-Anemic", `Iron status` = "Not sufficient", `B12 status` = "Sufficient")     %>% bind_cols(lead_row_ci(subgroup(FALSE, TRUE,  "Sufficient"))),
  tibble(`Hb status` = "Non-Anemic", `Iron status` = "Not sufficient", `B12 status` = "Not sufficient") %>% bind_cols(lead_row_ci(subgroup(FALSE, TRUE,  "Not sufficient"))),
  tibble(`Hb status` = "Total Non-Anemic", `Iron status` = "", `B12 status` = "") %>% bind_cols(lead_row_ci(d144 %>% filter(anemic == FALSE)))
)
print(supp_table3)






# Figure 3: Spearman correlation between hemoglobin and each biomarker, by age

biomarker_long <- long_data %>%
  select(pid, age_months, hb, bodyiron, b12, lead) %>%
  pivot_longer(c(bodyiron, b12, lead), names_to = "biomarker", values_to = "value") %>%
  filter(!is.na(hb), !is.na(value)) %>%
  mutate(biomarker = factor(biomarker, levels = c("bodyiron","b12","lead"),
                            labels = c("Body Iron (mg/kg)","Vitamin B12 (pg/ml)","Lead (\u00b5g/dL)")),
         age_label = factor(age_months, levels = biomarker_ages, labels = paste0(biomarker_ages, "m")))

cor_results <- biomarker_long %>%
  group_by(age_label, biomarker) %>%
  filter(n() >= 3) %>%
  summarise(n = n(), rho = suppressWarnings(cor.test(hb, value, method = "spearman"))$estimate,
            p = suppressWarnings(cor.test(hb, value, method = "spearman"))$p.value, .groups = "drop") %>%
  mutate(p_adj = p.adjust(p, method = "BH"))
print(cor_results)

star <- function(p) if (is.na(p)) "" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else ""

fig3 <- ggplot(cor_results, aes(age_label, biomarker, fill = rho)) +
  geom_tile(colour = "grey30") +
  geom_text(aes(label = paste0(sprintf("%.2f", rho), sapply(p_adj, star))), fontface = "bold", size = 3) +
  scale_fill_gradient2(low = "#B2182B", mid = "white", high = "#2166AC", midpoint = 0, limits = c(-1, 1)) +
  labs(x = "Age", y = NULL, fill = "Spearman rho") +
  theme_minimal()
print(fig3)







# Supplementary Table 1.a: data completeness by visit
# Numerator for each biomarker is restricted to children who ALSO have a concurrent hemoglobin value at that visit, matching long_data's denominator convention (Table 2, Supp Tables 2-4).

hb_col   <- c(`7`="Hb7", `15`="Hb15", `24`="Hb24", `36`="Hb36", `60`="Hb60", `84`="Hb84", `108`="Hb108", `144`="Hb144")
assessed <- sapply(hb_col, function(x) sum(!is.na(data[[x]])))

var_map <- tribble(
  ~age, ~Variable, ~column,
  7,   "Ferritin", "Ferritin7",   7,   "TfR", "tfr7",         7,   "Body Iron", "bodyiron7",
  15,  "Ferritin", "Ferritin15",  15,  "TfR", "tfr15",        15,  "Body Iron", "bodyiron15",  15,  "Lead", "Lead_15",
  24,  "Ferritin", "Ferritin24",  24,  "TfR", "tfr24",        24,  "Body Iron", "bodyiron24",  24,  "Lead", "Lead_24",
  108, "Ferritin", "Ferritin108", 108, "TfR", "tfr108",       108, "Body Iron", "bodyiron108", 108, "Vitamin B12", "B12_9y",
  144, "Ferritin", "Ferritin144", 144, "TfR", "tfr144",       144, "Body Iron", "bodyiron144", 144, "Vitamin B12", "B12_12y", 144, "Lead", "Lead_12y"
) %>% filter(column %in% names(data))

biomarker_completeness <- var_map %>%
  rowwise() %>%
  mutate(hb_col_name = hb_col[as.character(age)],
         n = sum(!is.na(data[[column]]) & !is.na(data[[hb_col_name]])),
         assessed_n = assessed[as.character(age)],
         cell = sprintf("%d (%.1f%%)", n, 100 * n / assessed_n)) %>%
  ungroup() %>%
  mutate(Variable = factor(Variable, levels = c("Ferritin","TfR","Body Iron","Lead","Vitamin B12"))) %>%
  select(age, Variable, cell) %>%
  pivot_wider(names_from = Variable, values_from = cell)

supp_table1a <- tibble(age = panel_ages) %>%
  left_join(biomarker_completeness, by = "age") %>%
  mutate(Hemoglobin = assessed[as.character(age)], .after = age) %>%
  mutate(across(-c(age, Hemoglobin), ~ replace_na(., "")))
print(supp_table1a)








# Supplementary Table 1.b: MAR sensitivity check (retained vs. dropped)

attrition_base <- data %>%
  filter(!is.na(Hb7)) %>%
  mutate(sex_group = case_when(CAFSEX == 1 ~ "Boys", CAFSEX == 2 ~ "Girls", TRUE ~ NA_character_),
         boys_flag = sex_group == "Boys",
         anemic_7m = mapply(hb_cutoff, Hb7, 7, sex_group))

safe_pct <- function(x) { n <- sum(!is.na(x)); if (n == 0) "-" else sprintf("%.1f%% (%d/%d)", 100 * sum(x, na.rm = TRUE) / n, sum(x, na.rm = TRUE), n) }

supp_table1b <- lapply(list(c(15,"Hb15"), c(24,"Hb24"), c(36,"Hb36"), c(60,"Hb60"),
                            c(84,"Hb84"), c(108,"Hb108"), c(144,"Hb144")), function(x) {
                              age <- as.numeric(x[1]); col <- x[2]
                              d   <- attrition_base %>% mutate(has_data = !is.na(.data[[col]]))
                              ret <- d %>% filter(has_data); drp <- d %>% filter(!has_data)
                              tibble(
                                `Follow-up (months)` = age, `N Retained` = nrow(ret), `N Dropped` = nrow(drp),
                                `Boys % (Retained)` = safe_pct(ret$boys_flag), `Boys % (Dropped)` = safe_pct(drp$boys_flag),
                                `Sex P` = attrition_test(d$has_data, d$boys_flag),
                                `Mean Hb7 (Retained)` = safe_meansd(ret$Hb7), `Mean Hb7 (Dropped)` = safe_meansd(drp$Hb7),
                                `Hb7 P` = safe_ttest_p(ret$Hb7, drp$Hb7),
                                `Anemic % 7m (Retained)` = safe_pct(ret$anemic_7m), `Anemic % 7m (Dropped)` = safe_pct(drp$anemic_7m),
                                `Anemia P` = attrition_test(d$has_data, d$anemic_7m)
                              )
                            }) %>% bind_rows()
print(supp_table1b)








# Supplementary Table 4: height-for-age z-score (HAZ) by anemia status

haz_long <- data %>%
  select(pid, HAZ_7, HAZ_15, HAZ_24, HAZ_108, HAZ_144) %>%
  pivot_longer(-pid, names_to = "col", values_to = "haz") %>%
  mutate(age_months = as.numeric(sub("HAZ_", "", col)), haz = as.numeric(haz)) %>%
  filter(!is.na(haz)) %>%
  select(pid, age_months, haz)

haz_anemia <- long_data %>% select(pid, age_months, hb_status) %>% inner_join(haz_long, by = c("pid","age_months"))

supp_table4 <- lapply(c(7, 15, 24, 108, 144), function(age) {
  d <- haz_anemia %>% filter(age_months == age)
  a  <- d$haz[d$hb_status == "Anemic"]
  n_ <- d$haz[d$hb_status == "Non-Anemic"]
  bind_rows(
    tibble(`Age (Months)` = age, `Hb Status` = "Anemic",     N = length(a),  `Mean HAZ (SD)` = safe_meansd(a),  P = ""),
    tibble(`Age (Months)` = age, `Hb Status` = "Non-Anemic", N = length(n_), `Mean HAZ (SD)` = safe_meansd(n_),
           P = if (length(a) > 1 && length(n_) > 1) safe_ttest_p(a, n_) else "-")
  )
}) %>% bind_rows()
print(supp_table4)