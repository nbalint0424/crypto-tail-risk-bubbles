# ============================================================
#  01 — Data: daily prices, log-returns, descriptive statistics
# ============================================================
#  Source : CoinMarketCap daily closes via the {crypto2} package
#  Coins  : BTC, ETH, BNB, XRP, SOL (fixed CMC ids, see R/utils.R)
#
#  Outputs
#    data/prices_wide.rds, data/log_returns.rds
#    output/tables/descriptive_stats.csv
#    output/figures/prices.png, returns.png, distributions.png
# ============================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(crypto2); library(dplyr); library(tidyr)
  library(ggplot2); library(moments); library(scales)
})

# ── Download ─────────────────────────────────────────────────
coin_list <- crypto_list(only_active = TRUE) |> filter(id %in% CMC_IDS)

raw <- crypto_history(coin_list  = coin_list,
                      start_date = START_DATE,
                      end_date   = END_DATE,
                      interval   = "daily",
                      sleep      = 0,       # raise to 1 if rate-limited
                      single_id  = TRUE)    # required since crypto2 v2.0
message("Downloaded ", nrow(raw), " rows")

# ── Wide price panel (one column per coin) ───────────────────
prices_wide <- raw |>
  mutate(date = as.Date(timestamp)) |>
  select(date, symbol, close) |>
  distinct(date, symbol, .keep_all = TRUE) |>
  pivot_wider(names_from = symbol, values_from = close) |>
  arrange(date) |>
  filter(!is.na(BTC)) |>                      # anchor on BTC (longest history)
  select(date, all_of(COINS))

# r_t = log(P_t) - log(P_{t-1})
log_returns <- prices_wide |>
  mutate(across(-date, ~ c(NA_real_, diff(log(.))))) |>
  slice(-1)

saveRDS(prices_wide, data_path("prices_wide.rds"))
saveRDS(log_returns, data_path("log_returns.rds"))

# ── Descriptive statistics ───────────────────────────────────
desc_stats <- log_returns |>
  select(all_of(COINS)) |>
  summarise(across(everything(), list(
    obs      = ~ sum(!is.na(.)),
    mean     = ~ mean(., na.rm = TRUE),
    sd       = ~ sd(., na.rm = TRUE),
    min      = ~ min(., na.rm = TRUE),
    max      = ~ max(., na.rm = TRUE),
    skewness = ~ skewness(., na.rm = TRUE),
    ex_kurt  = ~ kurtosis(., na.rm = TRUE) - 3   # moments::kurtosis is raw kurtosis
  ), .names = "{.col}__{.fn}")) |>
  pivot_longer(everything(), names_to = c("coin", "stat"), names_sep = "__") |>
  pivot_wider(names_from = stat, values_from = value)

print(desc_stats, n = Inf)
write.csv(desc_stats, table_path("descriptive_stats.csv"), row.names = FALSE)

# ── Figures ──────────────────────────────────────────────────
prices_indexed <- prices_wide |>
  pivot_longer(-date, names_to = "coin", values_to = "price") |>
  group_by(coin) |>
  mutate(index = price / first(na.omit(price)) * 100) |>
  ungroup()

p_prices <- ggplot(prices_indexed, aes(date, index, color = coin)) +
  geom_line(linewidth = 0.55, alpha = 0.9) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = COIN_COLORS) +
  labs(title = "Price index (rebased to 100, log scale)",
       x = NULL, y = "Index (log scale)", color = NULL) +
  theme_minimal(base_size = 13) + theme(legend.position = "top")
ggsave(fig_path("prices.png"), p_prices, width = 12, height = 5, dpi = 150)

ret_long <- log_returns |>
  pivot_longer(-date, names_to = "coin", values_to = "ret") |>
  filter(!is.na(ret))

p_returns <- ggplot(ret_long, aes(date, ret, color = coin)) +
  geom_line(linewidth = 0.35, alpha = 0.8) +
  scale_color_manual(values = COIN_COLORS) +
  facet_wrap(~coin, ncol = 1, scales = "free_y") +
  labs(title = "Daily log-returns", x = NULL, y = "Log-return") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")
ggsave(fig_path("returns.png"), p_returns, width = 12, height = 14, dpi = 150)

# Each coin is compared with a Normal using its OWN mean and sd.
normal_overlay <- ret_long |>
  group_by(coin) |>
  summarise(mu = mean(ret), s = sd(ret),
            lo = quantile(ret, 0.001), hi = quantile(ret, 0.999), .groups = "drop") |>
  rowwise() |>
  reframe(coin = coin, x = seq(lo, hi, length.out = 300), dens = dnorm(x, mu, s))

p_dist <- ggplot(ret_long, aes(ret)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, fill = "steelblue", alpha = 0.6) +
  geom_line(data = normal_overlay, aes(x, dens), color = "red", linetype = "dashed") +
  facet_wrap(~coin, scales = "free", ncol = 2) +
  labs(title = "Return distributions vs fitted Normal (dashed)",
       x = "Log-return", y = "Density") +
  theme_minimal(base_size = 11)
ggsave(fig_path("distributions.png"), p_dist, width = 11, height = 9, dpi = 150)

message("01 complete")
