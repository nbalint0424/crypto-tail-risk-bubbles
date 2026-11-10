# ============================================================
#  02 — Bubble detection: Phillips, Shi & Yu (2015) GSADF / BSADF
# ============================================================
#  Test on LOG-PRICE levels (explosive root), not on returns:
#    Δy_t = α + β·y_{t-1} + Σ φ_k Δy_{t-k} + ε_t
#    H0: β = 0 (unit root)   H1: β > 0 (explosive)
#
#  GSADF     : sup ADF over all sub-windows  → "is there any bubble?"
#  BSADF_t   : backward sup ADF at each date → date-stamping
#  minw      : floor(T·(0.01 + 1.8/√T))       (PSY 2015)
#  Crit. val.: Monte Carlo, 2,000 reps, seed 123
#  Episode   : BSADF_t > 95% CV for at least max(3, floor(log T)) days;
#              only positive (upward-explosive) episodes are kept.
#
#  Outputs
#    data/radf.rds, data/mc_critical_values.rds, data/bubble_regimes.rds
#    output/tables/gsadf_results.csv, bubble_episodes.csv, regime_statistics.csv
#    output/figures/bubble_<COIN>.png, regime_moments_BTC.png
# ============================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(exuber); library(dplyr); library(tidyr)
  library(ggplot2); library(moments); library(patchwork)
})

prices_wide <- readRDS(data_path("prices_wide.rds"))
log_returns <- readRDS(data_path("log_returns.rds"))

log_prices <- prices_wide |> mutate(across(-date, log))

series_of <- function(coin) {
  log_prices |> filter(!is.na(.data[[coin]])) |> select(date, y = all_of(coin))
}

# ── Recursive ADF statistics ─────────────────────────────────
radf_list <- lapply(setNames(COINS, COINS), function(coin) {
  s <- series_of(coin)
  message(sprintf("%-4s T = %d  minw = %d", coin, nrow(s), psy_minw(nrow(s))))
  radf(s$y, minw = psy_minw(nrow(s)), lag = 1)
})
saveRDS(radf_list, data_path("radf.rds"))

# ── Monte Carlo critical values (one set per sample length, cached) ──
cv_file <- data_path("mc_critical_values.rds")
T_vec   <- sapply(COINS, function(coin) nrow(series_of(coin)))

mc_cv_list <- if (file.exists(cv_file)) readRDS(cv_file) else list()
for (T_n in unique(T_vec)) {
  key <- as.character(T_n)
  if (is.null(mc_cv_list[[key]])) {
    message("Simulating critical values for T = ", T_n, " (slow, cached afterwards)")
    mc_cv_list[[key]] <- radf_mc_cv(n = T_n, minw = psy_minw(T_n), nrep = 2000, seed = SEED)
  }
}
saveRDS(mc_cv_list, cv_file)
cv_for <- function(coin) mc_cv_list[[as.character(T_vec[[coin]])]]

# ── GSADF test table ─────────────────────────────────────────
gsadf_table <- bind_rows(lapply(COINS, function(coin) {
  cv   <- cv_for(coin)$gsadf_cv
  stat <- as.numeric(radf_list[[coin]]$gsadf)
  data.frame(coin  = coin, T = T_vec[[coin]],
             GSADF = round(stat, 4),
             CV_90 = round(cv[["90%"]], 4),
             CV_95 = round(cv[["95%"]], 4),
             CV_99 = round(cv[["99%"]], 4),
             sig   = ifelse(stat > cv[["99%"]], "***",
                     ifelse(stat > cv[["95%"]], "**",
                     ifelse(stat > cv[["90%"]], "*", ""))))
}))
print(gsadf_table)
write.csv(gsadf_table, table_path("gsadf_results.csv"), row.names = FALSE)

# ── Date-stamp episodes ──────────────────────────────────────
# NOTE: exuber::datestamp() reports Start / Peak / End as positions in the
# ORIGINAL input series (the minw + lag burn-in is already accounted for),
# so they index the date vector directly — no extra offset.
episodes <- bind_rows(lapply(COINS, function(coin) {
  min_dur <- max(3L, floor(log(T_vec[[coin]])))
  ds <- tidy(datestamp(radf_list[[coin]], cv = cv_for(coin), min_duration = min_dur))
  if (nrow(ds) == 0) return(NULL)
  dates <- series_of(coin)$date
  ds |>
    filter(Signal == "positive") |>
    transmute(coin  = coin,
              start = dates[Start], peak = dates[Peak], end = dates[End],
              duration_days = as.integer(end - start) + 1L)
}))
print(episodes, n = Inf)
write.csv(episodes, table_path("bubble_episodes.csv"), row.names = FALSE)

# ── Daily 0/1 regime indicator ───────────────────────────────
regime_df <- log_prices |> select(date)
for (coin in COINS) {
  flag <- integer(nrow(regime_df))
  ep   <- episodes |> filter(coin == !!coin)
  for (i in seq_len(nrow(ep))) {
    flag[regime_df$date >= ep$start[i] & regime_df$date <= ep$end[i]] <- 1L
  }
  regime_df[[paste0(coin, "_bubble")]] <- flag
}

bubble_regimes <- prices_wide |>
  left_join(log_returns |> rename_with(~ paste0(., "_ret"), -date), by = "date") |>
  left_join(regime_df, by = "date")
saveRDS(bubble_regimes, data_path("bubble_regimes.rds"))

message("Bubble days per coin:")
print(colSums(select(regime_df, -date)))

# ── Return statistics by regime ──────────────────────────────
ret_regime <- bubble_regimes |>
  select(date, ends_with("_ret"), ends_with("_bubble")) |>
  pivot_longer(-date, names_to = c("coin", ".value"), names_sep = "_") |>
  filter(!is.na(ret)) |>
  mutate(regime = factor(ifelse(bubble == 1L, "Bubble", "Normal"), c("Normal", "Bubble")))

regime_stats <- ret_regime |>
  group_by(coin, regime) |>
  summarise(n        = n(),
            mean_pct = mean(ret) * 100,
            sd_pct   = sd(ret) * 100,
            skewness = skewness(ret),
            ex_kurt  = kurtosis(ret) - 3,
            min_pct  = min(ret) * 100,
            max_pct  = max(ret) * 100,
            .groups  = "drop") |>
  mutate(across(where(is.numeric) & !n, ~ round(., 4)))
print(regime_stats, n = Inf)
write.csv(regime_stats, table_path("regime_statistics.csv"), row.names = FALSE)

# ── Figure: log price + BSADF vs 95% CV (PSY 2015, Fig. 1 style) ──
for (coin in COINS) {
  s      <- series_of(coin)
  bsadf  <- as.numeric(radf_list[[coin]]$bsadf[, 1])
  b_df   <- data.frame(date = tail(s$date, length(bsadf)), bsadf = bsadf)
  cv95   <- cv_for(coin)$gsadf_cv[["95%"]]
  shade  <- make_shade(regime_df, coin)

  p_top <- ggplot(s, aes(date, y)) + shade_layer(shade, 0.25) +
    geom_line(color = COIN_COLORS[[coin]], linewidth = 0.6) +
    coord_cartesian(xlim = range(s$date)) +
    labs(title = paste0(coin, ": log price with PSY bubble episodes shaded"),
         x = NULL, y = "Log price (USD)") +
    theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))

  p_bot <- ggplot(b_df, aes(date, bsadf)) + shade_layer(shade, 0.25) +
    geom_line(color = "steelblue", linewidth = 0.5) +
    geom_hline(yintercept = cv95, color = "red", linetype = "dashed") +
    annotate("text", x = min(b_df$date) + 60, y = cv95 * 1.08,
             label = "95% critical value", color = "red", size = 3.5, hjust = 0) +
    coord_cartesian(xlim = range(s$date)) +
    labs(x = NULL, y = "BSADF") + theme_minimal(base_size = 12)

  ggsave(fig_path(paste0("bubble_", coin, ".png")),
         p_top / p_bot + plot_layout(heights = c(2, 1)), width = 12, height = 7, dpi = 150)
}

# ── Figure: 60-day rolling moments, BTC ──────────────────────
roll_stat <- function(x, window, FUN) {
  zoo::rollapplyr(x, window, function(w) {
    w <- w[!is.na(w)]
    if (length(w) >= 0.8 * window) FUN(w) else NA_real_
  }, fill = NA)
}

btc_mom <- log_returns |>
  filter(!is.na(BTC)) |>
  transmute(date,
            skew = roll_stat(BTC, 60, skewness),
            kurt = roll_stat(BTC, 60, function(x) kurtosis(x) - 3)) |>
  filter(!is.na(skew))

shade_btc <- make_shade(regime_df, "BTC")
mom_panel <- function(var, lab, col) {
  ggplot(btc_mom, aes(date, .data[[var]])) + shade_layer(shade_btc, 0.25) +
    geom_line(color = col, linewidth = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(x = NULL, y = lab) + theme_minimal(base_size = 12)
}
ggsave(fig_path("regime_moments_BTC.png"),
       (mom_panel("skew", "Rolling skewness (60d)", "steelblue") /
        mom_panel("kurt", "Rolling excess kurtosis (60d)", "darkgreen")) +
         plot_annotation(title = "BTC: 60-day rolling moments, bubble episodes shaded"),
       width = 12, height = 7, dpi = 150)

message("02 complete")
