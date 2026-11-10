# ============================================================
#  05 — Machine-learning quantile models: Lasso QR and XGBoost
# ============================================================
#  Both models predict conditional quantiles of the NEXT day's return
#  directly, with no distributional assumption.
#
#    Lasso quantile regression : rqPen::rq.pen.cv (check-loss, L1 penalty,
#                                λ chosen per quantile by 5-fold CV)
#    XGBoost                   : objective "reg:quantileerror" (pinball loss)
#
#  Features at day t (target = r_{t+1}):
#    r_{t-1}, r_{t-2}, r_{t-5}, |r_t|, |r_{t-1}|,
#    5- and 22-day realised variance, 22-day skewness,
#    GARCH(1,1)-t one-step-ahead σ for day t (from 04), bubble flag at t-1
#
#  Rolling design: 500-day window, models refit every 100 days.
#  Quantile grid τ ∈ {0.001, 0.005, 0.01, 0.025, 0.05}
#    VaR_α = q_α
#    ES_α  ≈ (1/α) ∫_0^α q_u du, integrating the piecewise-linear
#            interpolation of the quantile grid (q held flat below 0.001).
#
#  Outputs
#    data/ml_features.rds, data/ml_forecasts.rds
#    output/figures/ml_var_BTC.png, garch_vs_xgb_BTC.png, xgb_importance_BTC.png
# ============================================================

source("R/utils.R")
suppressPackageStartupMessages({
  library(rqPen); library(xgboost); library(dplyr); library(tidyr)
  library(ggplot2); library(moments)
})

log_returns     <- readRDS(data_path("log_returns.rds"))
bubble_regimes  <- readRDS(data_path("bubble_regimes.rds"))
garch_forecasts <- readRDS(data_path("garch_forecasts.rds"))

REFIT_EVERY <- 100L
Q_LEVELS    <- c(0.001, 0.005, 0.010, 0.025, 0.050)
FEATURES    <- c("r_lag1", "r_lag2", "r_lag5", "rv5", "rv22", "skew22",
                 "abs_r", "abs_lag1", "garch_sigma", "bubble_lag")
XGB_PARAMS  <- list(booster = "gbtree", objective = "reg:quantileerror",
                    eta = 0.05, max_depth = 4, subsample = 0.8,
                    colsample_bytree = 0.8, min_child_weight = 5)

# ── Feature engineering (row t uses information up to t only) ──
build_features <- function(coin) {
  garch_sig <- garch_forecasts[[coin]][["sGARCH_std"]] |> select(date, garch_sigma = sigma)
  log_returns |>
    select(date, r = all_of(coin)) |>
    filter(!is.na(r)) |>
    mutate(r_lag1   = lag(r, 1),
           r_lag2   = lag(r, 2),
           r_lag5   = lag(r, 5),
           rv5      = zoo::rollapplyr(r, 5,  var,      fill = NA),
           rv22     = zoo::rollapplyr(r, 22, var,      fill = NA),
           skew22   = zoo::rollapplyr(r, 22, skewness, fill = NA),
           abs_r    = abs(r),
           abs_lag1 = abs(lag(r, 1)),
           r_fwd    = lead(r, 1),
           date_fwd = lead(date, 1)) |>
    left_join(garch_sig, by = "date") |>
    left_join(select(bubble_regimes, date, bubble = all_of(paste0(coin, "_bubble"))), by = "date") |>
    mutate(bubble_lag = lag(bubble, 1)) |>
    filter(!is.na(r_fwd), !is.na(r_lag5), !is.na(rv22), !is.na(skew22)) |>
    # GARCH σ only exists after its own 500-day burn-in; earlier rows get 0.
    replace_na(list(bubble_lag = 0L, garch_sigma = 0))
}

features_list <- lapply(setNames(COINS, COINS), build_features)
saveRDS(features_list, data_path("ml_features.rds"))

# ── Quantile grid → VaR and ES ───────────────────────────────
# Weights w such that ES_α ≈ Σ w_k q(u_k) for the grid points u_k ≤ α.
es_weights <- function(grid, alpha) {
  u <- grid[grid <= alpha]
  w <- numeric(length(u))
  w[1] <- u[1]                                  # flat tail below the first grid point
  for (k in seq_len(length(u) - 1)) {           # trapezoid between grid points
    h <- u[k + 1] - u[k]
    w[k]     <- w[k]     + h / 2
    w[k + 1] <- w[k + 1] + h / 2
  }
  w / alpha
}

risk_from_quantiles <- function(q_mat) {
  out <- list()
  for (a in ALPHAS) {
    suffix <- sprintf("%02d", round(a * 100))
    idx    <- which(Q_LEVELS <= a + 1e-12)
    out[[paste0("VaR_", suffix)]] <- q_mat[, which.min(abs(Q_LEVELS - a))]
    out[[paste0("ES_",  suffix)]] <- as.numeric(q_mat[, idx, drop = FALSE] %*% es_weights(Q_LEVELS, a))
  }
  as.data.frame(out)
}

# ── Model wrappers ───────────────────────────────────────────
fit_lasso <- function(X, y) {
  # rqPen standardises predictors, so columns that are constant within the
  # window (e.g. the bubble flag when no bubble occurs in it) are dropped.
  keep <- apply(X, 2, sd) > 0
  fit <- withCallingHandlers(
    rq.pen.cv(X[, keep, drop = FALSE], y, tau = Q_LEVELS, penalty = "LASSO", nfolds = 5),
    warning = function(w) {   # non-unique LP solutions are routine for quantile regression
      if (grepl("nonunique", conditionMessage(w))) invokeRestart("muffleWarning")
    })
  list(fit = fit, keep = keep)
}
predict_lasso <- function(model, X) {
  # sort = TRUE rearranges crossing quantiles (Chernozhukov et al. 2010)
  suppressWarnings(predict(model$fit, newx = X[, model$keep, drop = FALSE], sort = TRUE))
}

fit_xgb <- function(X, y) {
  d <- xgb.DMatrix(X, label = y)
  lapply(Q_LEVELS, function(tau) {
    xgb.train(params = c(XGB_PARAMS, quantile_alpha = tau), data = d, nrounds = 100, verbose = 0)
  })
}
predict_xgb <- function(models, X) {
  sort(vapply(models, function(m) as.numeric(predict(m, xgb.DMatrix(X))), numeric(1)))
}

# ── Rolling forecasts ────────────────────────────────────────
fc_file <- data_path("ml_forecasts.rds")
if (file.exists(fc_file)) {
  ml_forecasts <- readRDS(fc_file)
} else {
  set.seed(SEED)   # CV fold assignment in rq.pen.cv
  ml_forecasts <- lapply(setNames(COINS, COINS), function(coin) {
    df    <- features_list[[coin]]
    X_all <- as.matrix(df[, FEATURES]); storage.mode(X_all) <- "double"
    y_all <- df$r_fwd
    n_fc  <- nrow(df) - WINDOW_SIZE

    q_lasso <- q_xgb <- matrix(NA_real_, n_fc, length(Q_LEVELS))
    lasso_m <- xgb_m <- NULL
    n_fail  <- c(lasso = 0L, xgb = 0L)

    for (i in seq_len(n_fc)) {
      if ((i - 1) %% REFIT_EVERY == 0) {
        tr <- i:(i + WINDOW_SIZE - 1)
        lasso_m <- tryCatch(fit_lasso(X_all[tr, ], y_all[tr]),
                            error = function(e) { n_fail["lasso"] <<- n_fail["lasso"] + 1L; NULL })
        xgb_m   <- tryCatch(fit_xgb(X_all[tr, ], y_all[tr]),
                            error = function(e) { n_fail["xgb"] <<- n_fail["xgb"] + 1L; NULL })
        if (i == 1 || i %% 500 == 1) message(sprintf("%s  step %d / %d", coin, i, n_fc))
      }
      x_new <- X_all[i + WINDOW_SIZE, , drop = FALSE]
      if (!is.null(lasso_m)) q_lasso[i, ] <- predict_lasso(lasso_m, x_new)
      if (!is.null(xgb_m))   q_xgb[i, ]   <- predict_xgb(xgb_m, x_new)
    }
    if (any(n_fail > 0)) warning(coin, ": failed refits ", paste(names(n_fail), n_fail, collapse = ", "))

    # Forecasts are dated by the day of the return they target (t+1),
    # the same convention as the GARCH forecasts, so the two align by date.
    base <- data.frame(date = df$date_fwd[WINDOW_SIZE + seq_len(n_fc)],
                       ret  = y_all[WINDOW_SIZE + seq_len(n_fc)])
    list(lasso = cbind(base, risk_from_quantiles(q_lasso)),
         xgb   = cbind(base, risk_from_quantiles(q_xgb)))
  })
  saveRDS(ml_forecasts, fc_file)
}

# ── Figures (BTC) ────────────────────────────────────────────
btc_ml <- bind_rows(lapply(c("lasso", "xgb"), function(m) mutate(ml_forecasts$BTC[[m]], model = m)))
shade  <- make_shade(bubble_regimes, "BTC", from = min(btc_ml$date))

p_ml <- ggplot(btc_ml, aes(date)) + shade_layer(shade, 0.15) +
  geom_line(aes(y = ret * 100), color = "grey40", linewidth = 0.3, alpha = 0.7) +
  geom_line(aes(y = VaR_05 * 100, color = model, linetype = "5%"), linewidth = 0.6) +
  geom_line(aes(y = VaR_01 * 100, color = model, linetype = "1%"), linewidth = 0.6) +
  scale_color_manual(values = c(lasso = "#9467bd", xgb = "#e377c2"),
                     labels = c(lasso = "Lasso QR", xgb = "XGBoost")) +
  scale_linetype_manual(values = c("5%" = "solid", "1%" = "dashed")) +
  labs(title = "BTC: rolling ML VaR forecasts",
       x = NULL, y = "Daily log-return (%)", color = "Model", linetype = "VaR level") +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("ml_var_BTC.png"), p_ml, width = 12, height = 6, dpi = 150)

cmp <- bind_rows(
  garch_forecasts$BTC$gjrGARCH_std |> select(date, ret, VaR_05) |> mutate(model = "GJR-GARCH(t)"),
  ml_forecasts$BTC$xgb            |> select(date, ret, VaR_05) |> mutate(model = "XGBoost")) |>
  group_by(date) |> filter(n() == 2) |> ungroup()

p_cmp <- ggplot(cmp, aes(date)) +
  shade_layer(make_shade(bubble_regimes, "BTC", from = min(cmp$date)), 0.15) +
  geom_line(data = filter(cmp, model == "XGBoost"), aes(y = ret * 100),
            color = "grey40", linewidth = 0.3, alpha = 0.6) +
  geom_line(aes(y = VaR_05 * 100, color = model), linewidth = 0.6) +
  scale_color_manual(values = c("GJR-GARCH(t)" = "#d62728", XGBoost = "#e377c2")) +
  labs(title = "BTC: GJR-GARCH(t) vs XGBoost, 5% VaR",
       x = NULL, y = "Daily log-return (%)", color = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("garch_vs_xgb_BTC.png"), p_cmp, width = 12, height = 5, dpi = 150)

# Feature importance from a full-sample fit (descriptive only, not used for forecasting).
df_btc <- features_list$BTC
d_btc  <- xgb.DMatrix(as.matrix(df_btc[, FEATURES]), label = df_btc$r_fwd)
imp <- bind_rows(lapply(c(0.05, 0.01), function(tau) {
  m <- xgb.train(params = c(XGB_PARAMS, quantile_alpha = tau), data = d_btc, nrounds = 200, verbose = 0)
  as.data.frame(xgb.importance(model = m)) |> mutate(level = paste0("VaR ", tau * 100, "%"))
}))

p_imp <- imp |>
  group_by(level) |> slice_max(Gain, n = 8) |> ungroup() |>
  ggplot(aes(reorder(Feature, Gain), Gain, fill = level)) +
  geom_col(position = "dodge", alpha = 0.85) + coord_flip() +
  scale_fill_manual(values = c("VaR 5%" = "#e377c2", "VaR 1%" = "#9467bd")) +
  labs(title = "XGBoost feature importance (BTC, full sample)",
       subtitle = "Gain = reduction in pinball loss attributable to each feature",
       x = NULL, y = "Gain", fill = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave(fig_path("xgb_importance_BTC.png"), p_imp, width = 10, height = 6, dpi = 150)

message("05 complete")
