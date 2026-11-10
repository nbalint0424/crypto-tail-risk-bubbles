# ============================================================
#  Shared configuration and helper functions
#  Sourced at the top of every pipeline script.
# ============================================================

# ── Configuration ────────────────────────────────────────────
COINS       <- c("BTC", "ETH", "BNB", "XRP", "SOL")
CMC_IDS     <- c(BTC = 1L, ETH = 1027L, BNB = 1839L, XRP = 52L, SOL = 5426L)
START_DATE  <- "2016-01-01"
END_DATE    <- "2026-06-09"   # fixed so results are reproducible

ALPHAS      <- c(0.05, 0.01)  # VaR / ES tail probabilities
WINDOW_SIZE <- 500L           # rolling estimation window (days)
SEED        <- 123L

DATA_DIR    <- "data"            # intermediate .rds files (git-ignored)
TABLE_DIR   <- "output/tables"
FIG_DIR     <- "output/figures"
for (d in c(DATA_DIR, TABLE_DIR, FIG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

data_path  <- function(f) file.path(DATA_DIR,  f)
table_path <- function(f) file.path(TABLE_DIR, f)
fig_path   <- function(f) file.path(FIG_DIR,   f)

COIN_COLORS <- c(BTC = "#F7931A", ETH = "#627EEA", BNB = "#F0B90B",
                 XRP = "#346AA9", SOL = "#9945FF")
REGIME_COLORS <- c(Normal = "steelblue", Bubble = "tomato")


# ── Formatting ───────────────────────────────────────────────
sig_stars <- function(p) {
  dplyr::case_when(is.na(p) ~ "",
                   p < 0.01 ~ "***",
                   p < 0.05 ~ "**",
                   p < 0.10 ~ "*",
                   TRUE     ~ "")
}


# ── Regime helpers ───────────────────────────────────────────

# Convert a 0/1 indicator into a table of contiguous episodes.
# Works on the FULL daily series, so gaps between episodes are kept.
episodes_from_indicator <- function(dates, indicator) {
  indicator <- as.integer(indicator)
  indicator[is.na(indicator)] <- 0L
  if (!any(indicator == 1L)) {
    return(data.frame(start = as.Date(character()), end = as.Date(character()),
                      duration = integer()))
  }
  r      <- rle(indicator)
  ends   <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1) + 1L)
  keep   <- r$values == 1L
  data.frame(start    = dates[starts[keep]],
             end      = dates[ends[keep]],
             duration = r$lengths[keep])
}

# Shading rectangles for ggplot (one per bubble episode).
make_shade <- function(regimes, coin, from = NULL) {
  ep <- episodes_from_indicator(regimes$date, regimes[[paste0(coin, "_bubble")]])
  if (nrow(ep) == 0) return(NULL)
  shade <- data.frame(xmin = ep$start, xmax = ep$end, ymin = -Inf, ymax = Inf)
  if (!is.null(from)) shade <- shade[shade$xmax >= from, , drop = FALSE]
  if (nrow(shade) == 0) NULL else shade
}

shade_layer <- function(shade, alpha = 0.2, fill = "tomato") {
  if (is.null(shade)) return(NULL)
  ggplot2::geom_rect(data = shade,
                     ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                     inherit.aes = FALSE, fill = fill, alpha = alpha)
}


# ── Backtesting ──────────────────────────────────────────────
# All tests take a logical violation vector with NAs already removed.

# Kupiec (1995) proportion-of-failures test.  LR_uc ~ chi2(1).
kupiec_pof <- function(viol, alpha) {
  n <- length(viol); x <- sum(viol)
  if (n == 0 || x == 0 || x == n) return(list(stat = NA_real_, pval = NA_real_))
  p  <- x / n
  lr <- -2 * (x * log(alpha / p) + (n - x) * log((1 - alpha) / (1 - p)))
  list(stat = lr, pval = 1 - pchisq(lr, df = 1))
}

# Christoffersen (1998) independence test (first-order Markov).  LR_ind ~ chi2(1).
christoffersen_ind <- function(viol) {
  v <- as.integer(viol); n <- length(v)
  if (n < 4) return(list(stat = NA_real_, pval = NA_real_))
  prev <- v[-n]; curr <- v[-1]
  n00 <- sum(prev == 0 & curr == 0); n01 <- sum(prev == 0 & curr == 1)
  n10 <- sum(prev == 1 & curr == 0); n11 <- sum(prev == 1 & curr == 1)
  pi01 <- n01 / (n00 + n01)
  pi11 <- n11 / (n10 + n11)
  pi_  <- (n01 + n11) / (n - 1)
  if (any(!is.finite(c(pi01, pi11))) || any(c(pi_, pi01, pi11) %in% c(0, 1))) {
    return(list(stat = NA_real_, pval = NA_real_))
  }
  lr <- -2 * ((n00 + n10) * log(1 - pi_) + (n01 + n11) * log(pi_) -
              n00 * log(1 - pi01) - n01 * log(pi01) -
              n10 * log(1 - pi11) - n11 * log(pi11))
  list(stat = lr, pval = 1 - pchisq(lr, df = 1))
}

# Pinball (quantile) loss: mean of (r - q) * (tau - 1{r < q}).
pinball_loss <- function(ret, q_hat, tau) {
  ok <- !is.na(ret) & !is.na(q_hat)
  if (sum(ok) < 2) return(NA_real_)
  r <- ret[ok]; q <- q_hat[ok]
  mean((r - q) * (tau - (r < q)))
}

# Basel-style traffic light, adapted to a violation RATE.
# The regulatory Basel rule counts exceptions of a 99% VaR over 250 days
# (green 0-4, yellow 5-9, red 10+).  Here the bands are expressed as excess
# over the nominal rate: green <= alpha + 1pp, yellow <= alpha + 4pp, else red.
basel_style_light <- function(viol_rate, alpha) {
  if (is.na(viol_rate)) return(NA_character_)
  excess <- viol_rate - alpha
  if (excess <= 0.01) "Green" else if (excess <= 0.04) "Yellow" else "Red"
}
