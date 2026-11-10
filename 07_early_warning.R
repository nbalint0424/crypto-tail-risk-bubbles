# ============================================================
#  07 — Early warning: do VaR forecasts rise before a bubble starts?
# ============================================================
#  For every PSY-dated bubble episode (from 02) that has forecast coverage,
#  take the 90 days BEFORE the origination date and compare the 5% VaR
#  path of GJR-GARCH(t) and XGBoost (VaR expressed as a positive loss).
#
#    1. Mann–Kendall trend test on the 90-day pre-bubble VaR path
#    2. Level change: mean VaR over days [-30, -1] vs [-60, -31], in %
#    3. Lead time: first day in the window on which VaR exceeds the
#       75th percentile of its own 90-day pre-bubble path
#
#  Episodes that start within the first 500 days of a coin's history have
#  no out-of-sample forecasts and are excluded.
#
#  Outputs
#    output/tables/early_warning_trends.csv, early_warning_levels.csv,
#                  early_warning_leadtime.csv
#    output/figures/ew_leadup.png, ew_leadtime.png
# ============================================================

source("R/utils.R")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(Kendall) })

garch_forecasts <- readRDS(data_path("garch_forecasts.rds"))
ml_forecasts    <- readRDS(data_path("ml_forecasts.rds"))
episodes        <- read.csv(table_path("bubble_episodes.csv")) |>
  mutate(start = as.Date(start))

PRE_WINDOW <- 90L
MIN_DAYS   <- 30L

var_paths <- bind_rows(lapply(COINS, function(coin) {
  full_join(
    garch_forecasts[[coin]]$gjrGARCH_std |> transmute(date, garch = -VaR_05),
    ml_forecasts[[coin]]$xgb             |> transmute(date, xgb   = -VaR_05),
    by = "date") |>
    mutate(coin = coin)
}))

pre_bubble <- bind_rows(lapply(seq_len(nrow(episodes)), function(i) {
  ep <- episodes[i, ]
  d  <- var_paths |>
    filter(coin == ep$coin, date >= ep$start - PRE_WINDOW, date < ep$start) |>
    mutate(episode = i, orig_date = ep$start, day = as.integer(date - ep$start))
  if (sum(!is.na(d$garch) & !is.na(d$xgb)) < MIN_DAYS) NULL else d
}))

message(sprintf("%d of %d episodes have pre-bubble forecast coverage",
                n_distinct(pre_bubble$episode), nrow(episodes)))

# ── 1. Mann–Kendall trend ────────────────────────────────────
mk <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 5) return(c(tau = NA, p = NA))
  res <- MannKendall(x)
  c(tau = unname(res$tau), p = unname(res$sl))   # two-sided p-value
}

trends <- pre_bubble |>
  arrange(episode, date) |>
  group_by(coin, episode, orig_date) |>
  summarise(n_days    = n(),
            tau_garch = mk(garch)[["tau"]], p_garch = mk(garch)[["p"]],
            tau_xgb   = mk(xgb)[["tau"]],   p_xgb   = mk(xgb)[["p"]],
            .groups = "drop") |>
  mutate(sig_garch = sig_stars(p_garch), sig_xgb = sig_stars(p_xgb),
         rising_garch = tau_garch > 0 & p_garch < 0.10,
         rising_xgb   = tau_xgb   > 0 & p_xgb   < 0.10)
print(as.data.frame(trends))
write.csv(trends, table_path("early_warning_trends.csv"), row.names = FALSE)

# ── 2. Level change in the run-up ────────────────────────────
levels <- pre_bubble |>
  mutate(period = case_when(day >= -60 & day < -30 ~ "early",
                            day >= -30 & day <   0 ~ "late")) |>
  filter(!is.na(period)) |>
  group_by(coin, episode, orig_date, period) |>
  summarise(garch = mean(garch, na.rm = TRUE), xgb = mean(xgb, na.rm = TRUE), .groups = "drop") |>
  pivot_wider(names_from = period, values_from = c(garch, xgb)) |>
  mutate(delta_garch_pct = (garch_late / garch_early - 1) * 100,
         delta_xgb_pct   = (xgb_late   / xgb_early   - 1) * 100)
write.csv(levels, table_path("early_warning_levels.csv"), row.names = FALSE)

# ── 3. Lead time ─────────────────────────────────────────────
first_exceed <- function(x, day) {
  hit <- which(x > quantile(x, 0.75, na.rm = TRUE))
  if (length(hit) == 0) NA_integer_ else abs(day[hit[1]])
}
lead_times <- pre_bubble |>
  arrange(episode, date) |>
  group_by(coin, episode, orig_date) |>
  summarise(lead_garch = first_exceed(garch, day),
            lead_xgb   = first_exceed(xgb, day), .groups = "drop") |>
  mutate(xgb_earlier = lead_xgb > lead_garch)
write.csv(lead_times, table_path("early_warning_leadtime.csv"), row.names = FALSE)

message(sprintf("Significant upward trend (p < 0.10): GARCH %d, XGBoost %d of %d episodes",
                sum(trends$rising_garch, na.rm = TRUE), sum(trends$rising_xgb, na.rm = TRUE),
                nrow(trends)))

# ── Figures ──────────────────────────────────────────────────
MODEL_COLS <- c("GJR-GARCH(t)" = "#d62728", XGBoost = "#e377c2")

leadup <- pre_bubble |>
  group_by(coin, day) |>
  summarise(`GJR-GARCH(t)` = mean(garch, na.rm = TRUE),
            XGBoost        = mean(xgb,   na.rm = TRUE), .groups = "drop") |>
  pivot_longer(-c(coin, day), names_to = "model", values_to = "var")

p_lead <- ggplot(leadup, aes(day, var * 100, color = model)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_line(linewidth = 0.6, na.rm = TRUE) +
  scale_color_manual(values = MODEL_COLS) +
  scale_x_continuous(breaks = seq(-90, 0, 30), labels = function(x) paste0(x, "d")) +
  facet_wrap(~coin, scales = "free_y", ncol = 2) +
  labs(title = "5% VaR in the 90 days before bubble origination",
       subtitle = "Average over covered episodes per coin | day 0 = PSY origination date",
       x = "Days relative to origination", y = "VaR (loss, %)", color = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")
ggsave(fig_path("ew_leadup.png"), p_lead, width = 12, height = 9, dpi = 150)

p_time <- lead_times |>
  pivot_longer(c(lead_garch, lead_xgb), names_to = "model", values_to = "days") |>
  mutate(model = ifelse(model == "lead_garch", "GJR-GARCH(t)", "XGBoost"),
         label = paste0(coin, "\n", orig_date)) |>
  filter(!is.na(days)) |>
  ggplot(aes(label, days, fill = model)) +
  geom_col(position = "dodge", alpha = 0.85) +
  scale_fill_manual(values = MODEL_COLS) +
  labs(title = "Early-warning lead time by episode",
       subtitle = "Days before origination at which VaR first exceeds the 75th percentile of its pre-bubble path",
       x = NULL, y = "Lead time (days)", fill = NULL) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", axis.text.x = element_text(size = 7, angle = 45, hjust = 1))
ggsave(fig_path("ew_leadtime.png"), p_time, width = 13, height = 5, dpi = 150)

message("07 complete")
