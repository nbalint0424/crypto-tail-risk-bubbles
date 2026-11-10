# ============================================================
#  03 — Return distributions: bubble vs normal regime
# ============================================================
#  1. Descriptive statistics incl. empirical VaR / ES by regime
#  2. Two-sample Kolmogorov–Smirnov test (bubble vs normal)
#  3. Jarque–Bera normality test per regime
#  4. Normal vs Student-t fit (MASS::fitdistr, compared by AIC)
#  5. Figures: densities, QQ plots, left-tail ECDF, volatility vs bubbles
#
#  Outputs
#    output/tables/regime_dist_stats.csv, ks_test_results.csv,
#                  jb_test_results.csv, dist_fit_results.csv
#    output/figures/densities.png, qq_regime.png, left_tail.png, vol_vs_bubble_BTC.png
# ============================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(moments)
})
# MASS is used via MASS:: only, so it never masks dplyr::select().

bubble_regimes <- readRDS(data_path("bubble_regimes.rds"))

ret_regime <- bubble_regimes |>
  select(date, ends_with("_ret"), ends_with("_bubble")) |>
  pivot_longer(-date, names_to = c("coin", ".value"), names_sep = "_") |>
  filter(!is.na(ret), !is.na(bubble)) |>
  mutate(regime = factor(ifelse(bubble == 1L, "Bubble", "Normal"), c("Normal", "Bubble")))

print(table(ret_regime$coin, ret_regime$regime))

# ── 1. Statistics by regime ──────────────────────────────────
emp_es <- function(x, a) mean(x[x <= quantile(x, a)])

regime_dist_stats <- ret_regime |>
  group_by(coin, regime) |>
  summarise(n          = n(),
            mean_pct   = mean(ret) * 100,
            median_pct = median(ret) * 100,
            sd_pct     = sd(ret) * 100,
            skewness   = skewness(ret),
            ex_kurt    = kurtosis(ret) - 3,
            VaR_5pct   = quantile(ret, 0.05) * 100,
            VaR_1pct   = quantile(ret, 0.01) * 100,
            ES_5pct    = emp_es(ret, 0.05) * 100,
            ES_1pct    = emp_es(ret, 0.01) * 100,
            .groups    = "drop") |>
  mutate(across(where(is.double), ~ round(., 4)))
print(regime_dist_stats, n = Inf)
write.csv(regime_dist_stats, table_path("regime_dist_stats.csv"), row.names = FALSE)

# ── 2. Kolmogorov–Smirnov: H0 same distribution ──────────────
ks_results <- bind_rows(lapply(COINS, function(cn) {
  xn <- ret_regime$ret[ret_regime$coin == cn & ret_regime$regime == "Normal"]
  xb <- ret_regime$ret[ret_regime$coin == cn & ret_regime$regime == "Bubble"]
  if (length(xb) < 5) {
    return(data.frame(coin = cn, n_normal = length(xn), n_bubble = length(xb),
                      D = NA, p_value = NA))
  }
  ks <- suppressWarnings(ks.test(xn, xb))   # ties in daily returns trigger a warning
  data.frame(coin = cn, n_normal = length(xn), n_bubble = length(xb),
             D = round(unname(ks$statistic), 4), p_value = signif(ks$p.value, 4))
})) |> mutate(sig = sig_stars(p_value))
print(ks_results)
write.csv(ks_results, table_path("ks_test_results.csv"), row.names = FALSE)

# ── 3. Jarque–Bera per regime ────────────────────────────────
jb_results <- ret_regime |>
  group_by(coin, regime) |>
  filter(n() >= 10) |>
  summarise(n       = n(),
            JB_stat = round(unname(jarque.test(ret)$statistic), 2),
            p_value = signif(jarque.test(ret)$p.value, 4),
            .groups = "drop") |>
  mutate(sig = sig_stars(p_value))
print(jb_results, n = Inf)
write.csv(jb_results, table_path("jb_test_results.csv"), row.names = FALSE)

# ── 4. Normal vs Student-t (lower AIC wins; smaller df = fatter tails) ──
fit_one <- function(x) {
  fn <- tryCatch(MASS::fitdistr(x, "normal"), error = function(e) NULL)
  ft <- tryCatch(suppressWarnings(MASS::fitdistr(x, "t")), error = function(e) NULL)
  aic_n <- if (is.null(fn)) NA else AIC(fn)
  aic_t <- if (is.null(ft)) NA else AIC(ft)
  data.frame(AIC_N = round(aic_n, 2), AIC_t = round(aic_t, 2),
             df_t  = if (is.null(ft)) NA else round(unname(ft$estimate["df"]), 2),
             best_fit = ifelse(is.na(aic_n) | is.na(aic_t), NA,
                               ifelse(aic_t < aic_n, "Student-t", "Normal")))
}

fit_results <- ret_regime |>
  group_by(coin, regime) |>
  filter(n() >= 20) |>
  group_modify(~ cbind(n = nrow(.x), fit_one(.x$ret))) |>
  ungroup()
print(fit_results, n = Inf)
write.csv(fit_results, table_path("dist_fit_results.csv"), row.names = FALSE)

# ── 5. Figures ───────────────────────────────────────────────
p_dens <- ggplot(ret_regime, aes(ret * 100, fill = regime, color = regime)) +
  geom_density(alpha = 0.35, linewidth = 0.7) +
  scale_fill_manual(values = REGIME_COLORS) + scale_color_manual(values = REGIME_COLORS) +
  facet_wrap(~coin, scales = "free", ncol = 2) +
  coord_cartesian(xlim = c(-20, 20)) +
  labs(title = "Return density by regime", x = "Daily log-return (%)", y = "Density",
       fill = NULL, color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("densities.png"), p_dens, width = 11, height = 9, dpi = 150)

qq_data <- ret_regime |>
  group_by(coin, regime) |>
  filter(n() >= 5) |>
  mutate(theoretical = qnorm(ppoints(n())),
         sample_q    = sort((ret - mean(ret)) / sd(ret))) |>
  ungroup()

p_qq <- ggplot(qq_data, aes(theoretical, sample_q, color = regime)) +
  geom_point(size = 0.4, alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_manual(values = REGIME_COLORS) +
  facet_wrap(~coin, scales = "free", ncol = 2) +
  labs(title = "QQ plot vs Normal, by regime",
       x = "Theoretical Normal quantiles", y = "Standardised sample quantiles", color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("qq_regime.png"), p_qq, width = 11, height = 9, dpi = 150)

tail_data <- ret_regime |>
  group_by(coin, regime) |>
  arrange(ret, .by_group = TRUE) |>
  mutate(ecdf_p = seq_along(ret) / (n() + 1)) |>
  filter(ecdf_p <= 0.10) |>
  ungroup()

p_tail <- ggplot(tail_data, aes(ret * 100, ecdf_p, color = regime)) +
  geom_line(linewidth = 0.8) +
  scale_color_manual(values = REGIME_COLORS) +
  facet_wrap(~coin, scales = "free_x", ncol = 2) +
  labs(title = "Left tail: empirical CDF of the worst 10% of days",
       x = "Daily log-return (%)", y = "Empirical probability", color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("left_tail.png"), p_tail, width = 11, height = 9, dpi = 150)

# Sanity check: PSY flags explosive price growth, not simply high volatility.
btc_vol <- bubble_regimes |>
  filter(!is.na(BTC_ret)) |>
  transmute(date, vol = zoo::rollapplyr(BTC_ret, 30, sd, fill = NA) * sqrt(365) * 100) |>
  filter(!is.na(vol))

p_vol <- ggplot(btc_vol, aes(date, vol)) +
  shade_layer(make_shade(bubble_regimes, "BTC"), 0.25) +
  geom_line(color = COIN_COLORS[["BTC"]], linewidth = 0.6) +
  labs(title = "BTC: 30-day annualised volatility, bubble episodes shaded",
       subtitle = "Crashes are high-volatility periods without explosive price dynamics",
       x = NULL, y = "Annualised volatility (%)") +
  theme_minimal(base_size = 12)
ggsave(fig_path("vol_vs_bubble_BTC.png"), p_vol, width = 12, height = 5, dpi = 150)

message("03 complete")
