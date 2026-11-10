# ============================================================
#  04 — GARCH benchmarks: estimation and rolling VaR / ES forecasts
# ============================================================
#  Variance models (all (1,1), constant mean):
#    sGARCH   σ²_t = ω + α ε²_{t-1} + β σ²_{t-1}                        (Bollerslev 1986)
#    eGARCH   log σ²_t = ω + α|z_{t-1}| + γ z_{t-1} + β log σ²_{t-1}     (Nelson 1991)
#    gjrGARCH σ²_t = ω + (α + γ·1{ε_{t-1}<0}) ε²_{t-1} + β σ²_{t-1}      (GJR 1993)
#  Innovations: Normal (norm), Student-t (std), skewed Student-t (sstd)
#    → 9 specifications, all fitted full-sample for the AIC/BIC comparison.
#
#  Rolling forecasts (t and skew-t only, Normal is dominated on AIC):
#    500-day moving window, refit every 25 days, 1-day-ahead.
#    VaR_α = μ + σ · q_α(Z)
#    ES_α  = μ + σ · ES_α(Z),  ES_α(Z) = (1/α) ∫_0^α q_u(Z) du
#    where q(Z) uses the shape / skew parameters estimated in that window.
#
#  Outputs
#    data/garch_fits.rds, data/garch_forecasts.rds
#    output/tables/garch_ic_table.csv, garch_fit_table.csv
#    output/figures/garch_vol_BTC.png, garch_var_BTC.png, garch_vol_by_regime_BTC.png
# ============================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(rugarch); library(dplyr); library(tidyr); library(ggplot2)
})

log_returns    <- readRDS(data_path("log_returns.rds"))
bubble_regimes <- readRDS(data_path("bubble_regimes.rds"))

REFIT_EVERY <- 25L

# ── Specifications ───────────────────────────────────────────
make_spec <- function(model, dist) {
  ugarchspec(variance.model     = list(model = model, garchOrder = c(1, 1)),
             mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
             distribution.model = dist)
}
spec_grid   <- expand.grid(model = c("sGARCH", "eGARCH", "gjrGARCH"),
                           dist  = c("norm", "std", "sstd"), stringsAsFactors = FALSE)
garch_specs <- setNames(Map(make_spec, spec_grid$model, spec_grid$dist),
                        paste(spec_grid$model, spec_grid$dist, sep = "_"))
fc_names    <- grep("_(std|sstd)$", names(garch_specs), value = TRUE)

returns_of <- function(coin) {
  log_returns |> filter(!is.na(.data[[coin]])) |> select(date, r = all_of(coin))
}

# ── Full-sample estimation (for model comparison and plots) ──
fits_file <- data_path("garch_fits.rds")
if (file.exists(fits_file)) {
  garch_fits <- readRDS(fits_file)
} else {
  garch_fits <- lapply(setNames(COINS, COINS), function(coin) {
    r <- returns_of(coin)$r
    lapply(garch_specs, function(sp) {
      tryCatch(ugarchfit(sp, data = r, solver = "hybrid"),
               error = function(e) { warning(coin, ": ", conditionMessage(e)); NULL })
    })
  })
  saveRDS(garch_fits, fits_file)
}

ic_table <- bind_rows(lapply(COINS, function(coin) {
  bind_rows(lapply(names(garch_specs), function(sn) {
    fit <- garch_fits[[coin]][[sn]]
    if (is.null(fit) || fit@fit$convergence != 0) return(NULL)
    ic <- infocriteria(fit)
    data.frame(coin = coin, model = sn,
               AIC = round(ic[1], 4), BIC = round(ic[2], 4),
               LogLik = round(likelihood(fit), 2))
  }))
}))
write.csv(ic_table, table_path("garch_ic_table.csv"), row.names = FALSE)
print(ic_table |> group_by(coin) |> slice_min(AIC, n = 1, with_ties = FALSE))

param_table <- bind_rows(lapply(COINS, function(coin) {
  bind_rows(lapply(names(garch_specs), function(sn) {
    fit <- garch_fits[[coin]][[sn]]
    if (is.null(fit) || fit@fit$convergence != 0) return(NULL)
    data.frame(coin = coin, model = sn, t(round(coef(fit), 6)), check.names = FALSE)
  }))
}))
write.csv(param_table, table_path("garch_fit_table.csv"), row.names = FALSE)

# ── Innovation quantiles and expected shortfall ──────────────
# rugarch::qdist uses the standardised (zero-mean, unit-variance) parameterisation.
z_quantile <- function(p, dist, skew, shape) {
  qdist(distribution = dist, p = p, mu = 0, sigma = 1,
        skew  = ifelse(is.na(skew), 1, skew),
        shape = ifelse(is.na(shape), 5, shape))
}

# ES_α(Z) = (1/α) ∫_0^α q_u du, evaluated with a midpoint rule on (0, α).
z_es <- function(alpha, dist, skew, shape, n_grid = 2000L) {
  u <- alpha * (seq_len(n_grid) - 0.5) / n_grid
  mean(z_quantile(u, dist, skew, shape))
}

# Shape/skew only change at refits, so evaluate each parameter set once.
risk_from_density <- function(dens, dist) {
  key    <- paste(dens$Skew, dens$Shape)
  params <- unique(data.frame(key, skew = dens$Skew, shape = dens$Shape))
  out    <- data.frame(row.names = seq_len(nrow(dens)))
  for (a in ALPHAS) {
    suffix <- sprintf("%02d", round(a * 100))
    q_z  <- setNames(mapply(z_quantile, a, dist, params$skew, params$shape), params$key)
    es_z <- setNames(mapply(z_es,       a, dist, params$skew, params$shape), params$key)
    out[[paste0("VaR_", suffix)]] <- dens$Mu + dens$Sigma * q_z[key]
    out[[paste0("ES_",  suffix)]] <- dens$Mu + dens$Sigma * es_z[key]
  }
  out
}

# ── Rolling-window forecasts ─────────────────────────────────
fc_file <- data_path("garch_forecasts.rds")
if (file.exists(fc_file)) {
  garch_forecasts <- readRDS(fc_file)
} else {
  garch_forecasts <- lapply(setNames(COINS, COINS), function(coin) {
    rd   <- returns_of(coin)
    n_fc <- nrow(rd) - WINDOW_SIZE
    lapply(setNames(fc_names, fc_names), function(sn) {
      message(sprintf("%-4s %-14s %d forecasts", coin, sn, n_fc))
      roll <- tryCatch(
        ugarchroll(garch_specs[[sn]], data = rd$r, n.ahead = 1,
                   forecast.length = n_fc, refit.every = REFIT_EVERY,
                   refit.window = "moving", window.size = WINDOW_SIZE,
                   solver = "hybrid", calculate.VaR = FALSE, keep.coef = FALSE),
        error = function(e) { warning(coin, " ", sn, ": ", conditionMessage(e)); NULL })
      if (is.null(roll)) return(NULL)

      # Windows that failed to converge are re-estimated with a global solver.
      if (convergence(roll) != 0) roll <- resume(roll, solver = "gosolnp")
      if (convergence(roll) != 0) {
        warning(coin, " ", sn, ": some windows still did not converge; their forecasts are NA")
      }

      dens <- as.data.frame(roll@forecast$density)
      cbind(data.frame(date  = tail(rd$date, n_fc),
                       ret   = dens$Realized,
                       mu    = dens$Mu,
                       sigma = dens$Sigma),
            risk_from_density(dens, garch_specs[[sn]]@model$modeldesc$distribution))
    })
  })
  saveRDS(garch_forecasts, fc_file)
}

# ── Figures (BTC) ────────────────────────────────────────────
best_btc <- ic_table |> filter(coin == "BTC") |> slice_min(AIC, n = 1) |> pull(model)
fit_btc  <- garch_fits[["BTC"]][[best_btc]]
btc_sig  <- returns_of("BTC") |>
  mutate(vol = as.numeric(sigma(fit_btc)) * sqrt(365) * 100) |>
  left_join(select(bubble_regimes, date, bubble = BTC_bubble), by = "date")
shade_btc <- make_shade(bubble_regimes, "BTC")

p_vol <- ggplot(btc_sig, aes(date, vol)) + shade_layer(shade_btc) +
  geom_line(color = COIN_COLORS[["BTC"]], linewidth = 0.5) +
  labs(title = paste0("BTC: conditional volatility, ", best_btc, " (best AIC)"),
       subtitle = "Annualised (365 days); bubble episodes shaded",
       x = NULL, y = "Annualised volatility (%)") +
  theme_minimal(base_size = 12)
ggsave(fig_path("garch_vol_BTC.png"), p_vol, width = 12, height = 5, dpi = 150)

btc_fc <- bind_rows(lapply(grep("_std$", fc_names, value = TRUE), function(sn) {
  garch_forecasts[["BTC"]][[sn]] |> mutate(model = sub("_std$", "", sn))
}))

p_var <- ggplot(btc_fc, aes(date)) +
  shade_layer(make_shade(bubble_regimes, "BTC", from = min(btc_fc$date)), 0.15) +
  geom_line(aes(y = ret * 100), color = "grey40", linewidth = 0.3, alpha = 0.7) +
  geom_line(aes(y = VaR_05 * 100, color = model, linetype = "5%"), linewidth = 0.6) +
  geom_line(aes(y = VaR_01 * 100, color = model, linetype = "1%"), linewidth = 0.6) +
  scale_color_manual(values = c(sGARCH = "#1f77b4", eGARCH = "#2ca02c", gjrGARCH = "#d62728")) +
  scale_linetype_manual(values = c("5%" = "solid", "1%" = "dashed")) +
  labs(title = "BTC: rolling 1-day-ahead VaR, Student-t GARCH",
       subtitle = "Grey = realised returns | 500-day window, refit every 25 days | bubble episodes shaded",
       x = NULL, y = "Daily log-return (%)", color = "Model", linetype = "VaR level") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("garch_var_BTC.png"), p_var, width = 12, height = 6, dpi = 150)

p_box <- btc_sig |>
  filter(!is.na(bubble)) |>
  mutate(regime = factor(ifelse(bubble == 1L, "Bubble", "Normal"), c("Normal", "Bubble"))) |>
  ggplot(aes(regime, vol, fill = regime)) +
  geom_violin(alpha = 0.4) + geom_boxplot(width = 0.15, outlier.size = 0.5) +
  scale_fill_manual(values = REGIME_COLORS) +
  labs(title = paste0("BTC: conditional volatility by regime (", best_btc, ")"),
       x = NULL, y = "Annualised volatility (%)") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(fig_path("garch_vol_by_regime_BTC.png"), p_box, width = 8, height = 6, dpi = 150)

message("04 complete")
