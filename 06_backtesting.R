# ============================================================
#  06 — VaR backtesting, full sample and by regime
# ============================================================
#  For every coin × model × α ∈ {5%, 1%}:
#    Kupiec (1995) POF test               H0: violation rate = α          LR ~ χ²(1)
#    Christoffersen (1998) independence   H0: violations not clustered    LR ~ χ²(1)
#    Conditional coverage                 LR_cc = LR_uc + LR_ind          LR ~ χ²(2)
#    Basel-style traffic light            see basel_style_light() in R/utils.R
#    Pinball loss                         mean (r − q)(τ − 1{r < q}); lower is better
#  Rates, Kupiec and pinball loss are also computed separately on normal-regime
#  and bubble-regime days. Christoffersen is full-sample only, since splitting
#  the series breaks the day-to-day sequence the test relies on.
#
#  Outputs
#    output/tables/backtest_results.csv, backtest_summary.csv, backtest_regime_loss.csv
#    output/figures/violation_heatmap.png, pinball_by_regime.png, rolling_violations_BTC.png
# ============================================================

source("R/utils.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2) })

bubble_regimes  <- readRDS(data_path("bubble_regimes.rds"))
garch_forecasts <- readRDS(data_path("garch_forecasts.rds"))
ml_forecasts    <- readRDS(data_path("ml_forecasts.rds"))

MODEL_LABELS <- c(sGARCH_std = "sGARCH(t)", eGARCH_std = "eGARCH(t)", gjrGARCH_std = "GJR-GARCH(t)",
                  sGARCH_sstd = "sGARCH(skew-t)", eGARCH_sstd = "eGARCH(skew-t)",
                  gjrGARCH_sstd = "GJR-GARCH(skew-t)", xgb = "XGBoost", lasso = "Lasso QR")
SUMMARY_MODELS <- c("sGARCH(t)", "eGARCH(t)", "GJR-GARCH(t)", "XGBoost", "Lasso QR")

# ── One long table of forecasts with the regime of each day ──
regime_long <- bubble_regimes |>
  select(date, ends_with("_bubble")) |>
  pivot_longer(-date, names_to = "coin", names_pattern = "(.*)_bubble", values_to = "bubble")

all_fc <- bind_rows(
  lapply(COINS, function(coin) bind_rows(c(
    lapply(names(garch_forecasts[[coin]]), function(m)
      garch_forecasts[[coin]][[m]] |> mutate(coin = coin, model = m)),
    lapply(names(ml_forecasts[[coin]]), function(m)
      ml_forecasts[[coin]][[m]] |> mutate(coin = coin, model = m)))))) |>
  select(coin, model, date, ret, VaR_05, VaR_01) |>
  left_join(regime_long, by = c("coin", "date")) |>
  mutate(model  = unname(MODEL_LABELS[model]),
         bubble = coalesce(as.integer(bubble), 0L))

# ── Backtest one forecast series at one level ────────────────
backtest_one <- function(ret, q_hat, bubble, alpha) {
  ok <- !is.na(ret) & !is.na(q_hat)
  ret <- ret[ok]; q_hat <- q_hat[ok]; bubble <- bubble[ok]
  viol <- ret < q_hat

  by_regime <- function(sel, tag) {
    v  <- viol[sel]
    kp <- if (length(v) > 5) kupiec_pof(v, alpha) else list(pval = NA_real_)
    vr <- if (length(v) > 0) mean(v) else NA_real_
    setNames(data.frame(length(v), vr * 100, basel_style_light(vr, alpha),
                        kp$pval, pinball_loss(ret[sel], q_hat[sel], alpha)),
             paste0(c("n_", "viol_", "light_", "kupiec_p_", "pinball_"), tag))
  }

  kp  <- kupiec_pof(viol, alpha)
  ind <- christoffersen_ind(viol)
  cc  <- kp$stat + ind$stat
  cbind(
    data.frame(n_obs = length(viol), viol_rate = mean(viol) * 100,
               light = basel_style_light(mean(viol), alpha),
               kupiec_stat = kp$stat, kupiec_p = kp$pval,
               ind_stat = ind$stat, ind_p = ind$pval,
               cc_stat = cc, cc_p = 1 - pchisq(cc, df = 2),
               pinball = pinball_loss(ret, q_hat, alpha)),
    by_regime(bubble == 0L, "normal"),
    by_regime(bubble == 1L, "bubble"))
}

backtest_results <- all_fc |>
  group_by(coin, model) |>
  group_modify(~ bind_rows(
    cbind(alpha = 0.05, backtest_one(.x$ret, .x$VaR_05, .x$bubble, 0.05)),
    cbind(alpha = 0.01, backtest_one(.x$ret, .x$VaR_01, .x$bubble, 0.01)))) |>
  ungroup() |>
  mutate(across(where(is.double) & !alpha, ~ signif(., 5)))
write.csv(backtest_results, table_path("backtest_results.csv"), row.names = FALSE)

summary_table <- backtest_results |>
  filter(model %in% SUMMARY_MODELS) |>
  transmute(coin, model, alpha, viol_rate = round(viol_rate, 2), light,
            kupiec = sig_stars(kupiec_p), indep = sig_stars(ind_p),
            viol_normal = round(viol_normal, 2), viol_bubble = round(viol_bubble, 2),
            n_bubble, pinball_normal, pinball_bubble) |>
  arrange(alpha, coin, model)
print(filter(summary_table, alpha == 0.05), n = Inf)
write.csv(summary_table, table_path("backtest_summary.csv"), row.names = FALSE)

regime_loss <- backtest_results |>
  filter(model %in% SUMMARY_MODELS) |>
  transmute(coin, model, alpha, pinball_normal, pinball_bubble,
            loss_diff = pinball_bubble - pinball_normal)
write.csv(regime_loss, table_path("backtest_regime_loss.csv"), row.names = FALSE)

# ── Figures ──────────────────────────────────────────────────
heat <- summary_table |>
  mutate(level = factor(paste0("VaR ", alpha * 100, "%"), c("VaR 5%", "VaR 1%")),
         model = factor(model, rev(SUMMARY_MODELS)), coin = factor(coin, COINS))

p_heat <- ggplot(heat, aes(coin, model, fill = light)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f%%", viol_rate)), size = 3.2) +
  scale_fill_manual(values = c(Green = "#9be3b0", Yellow = "#ffe08a", Red = "#f4a3a3"),
                    na.value = "grey90", name = "Basel-style\nlight") +
  facet_wrap(~level, ncol = 1) +
  labs(title = "Full-sample VaR violation rates",
       subtitle = "Green: rate <= alpha + 1pp | Yellow: <= alpha + 4pp | Red: above",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11)
ggsave(fig_path("violation_heatmap.png"), p_heat, width = 9, height = 8, dpi = 150)

p_loss <- regime_loss |>
  filter(alpha == 0.05) |>
  pivot_longer(c(pinball_normal, pinball_bubble), names_to = "regime", values_to = "loss",
               names_prefix = "pinball_") |>
  mutate(regime = factor(tools::toTitleCase(regime), c("Normal", "Bubble")),
         model  = factor(model, SUMMARY_MODELS)) |>
  filter(!is.na(loss)) |>
  ggplot(aes(model, loss, fill = regime)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = REGIME_COLORS) +
  facet_wrap(~coin, scales = "free_y", ncol = 2) +
  labs(title = "Pinball loss by regime, 5% VaR (lower is better)",
       x = NULL, y = "Average pinball loss", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(fig_path("pinball_by_regime.png"), p_loss, width = 12, height = 10, dpi = 150)

roll_viol <- all_fc |>
  filter(coin == "BTC", model %in% c("GJR-GARCH(t)", "XGBoost"), !is.na(VaR_05)) |>
  arrange(model, date) |>
  group_by(model) |>
  mutate(rate = zoo::rollmeanr(as.numeric(ret < VaR_05), 30, fill = NA) * 100) |>
  ungroup()

p_roll <- ggplot(roll_viol, aes(date, rate, color = model)) +
  shade_layer(make_shade(bubble_regimes, "BTC", from = min(roll_viol$date)), 0.15) +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  geom_hline(yintercept = 5, linetype = "dashed") +
  scale_color_manual(values = c("GJR-GARCH(t)" = "#d62728", XGBoost = "#e377c2")) +
  labs(title = "BTC: rolling 30-day 5% VaR violation rate",
       subtitle = "Dashed line = nominal 5% | bubble episodes shaded",
       x = NULL, y = "Violation rate (%)", color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("rolling_violations_BTC.png"), p_roll, width = 12, height = 5, dpi = 150)

message("06 complete")
