# Tail Risk and Bubble Regimes in Cryptocurrency Markets

Code for my MSc Quantitative Finance thesis at WU Vienna, *Quantitative Modeling of Tail Risk and Bubble Regimes in Cryptocurrency Markets*.

**Question:** do standard Value-at-Risk models break down during speculative bubbles, and can machine-learning quantile models do better?

The pipeline:

1. date-stamps explosive bubble episodes in BTC, ETH, BNB, XRP and SOL (daily, 2016–2026) using the Phillips–Shi–Yu (2015) GSADF/BSADF test;
2. forecasts 1-day VaR and Expected Shortfall out of sample with nine GARCH specifications and two ML quantile models (Lasso quantile regression, XGBoost);
3. backtests every model on the full sample and separately on bubble vs normal days.

## Methods

| Step | Script | What it does |
|---|---|---|
| Data | `01_data.R` | CoinMarketCap daily closes via `crypto2`, log-returns, descriptive statistics |
| Bubble detection | `02_bubble_detection.R` | GSADF test on log prices; BSADF date-stamping against 2,000-replication Monte Carlo critical values; daily bubble indicator |
| Distributions | `03_regime_distributions.R` | Moments, empirical VaR/ES, Kolmogorov–Smirnov, Jarque–Bera, Normal vs Student-t fit by regime |
| GARCH | `04_garch.R` | sGARCH / eGARCH / GJR-GARCH × Normal / Student-t / skew-t; 500-day rolling window, refit every 25 days; parametric VaR and ES from the fitted innovation distribution |
| ML quantile models | `05_ml_quantile.R` | Lasso QR (`rqPen`, CV-tuned per quantile) and XGBoost (pinball-loss objective) on lagged returns, realised variance, skewness, GARCH σ and the bubble flag; ES by integrating the predicted quantile curve |
| Backtesting | `06_backtesting.R` | Kupiec POF, Christoffersen independence and conditional coverage, Basel-style traffic light, pinball loss, all split by regime |
| Early warning | `07_early_warning.R` | Does VaR rise in the 90 days before a bubble starts? Mann–Kendall trend, level change, lead time |

The backtests are written from scratch in [`R/utils.R`](R/utils.R) rather than taken from a package.

## Results

*Headline findings from the thesis (daily data to 9 June 2026).*

**Bubbles are detected in every coin.** The GSADF test rejects the no-bubble null at the 1% level for BTC, ETH, BNB and XRP and at 5% for SOL. Date-stamping finds 27 episodes, most of them short (8–30 days). The longest are BTC Nov 2017–Jan 2018, ETH and XRP Apr–Jul 2017, and BNB Feb–May 2021.

| Coin | T | GSADF | 95% CV | 99% CV |
|---|---:|---:|---:|---:|
| BTC | 3,813 | 4.03 | 2.52 | 2.95 |
| ETH | 3,813 | 3.75 | 2.52 | 2.95 |
| BNB | 3,242 | 5.13 | 2.49 | 2.94 |
| XRP | 3,813 | 7.43 | 2.52 | 2.95 |
| SOL | 2,252 | 2.58 | 2.43 | 2.93 |

**Student-t GARCH is well calibrated over the full sample; XGBoost is not.** At the 5% level, Student-t GARCH violation rates are 4.7–5.7% for BTC, ETH, BNB and XRP, with the Kupiec test not rejecting at the 5% level. XGBoost under-predicts risk for every coin (6.8–7.8% violations, Kupiec p < 0.001). GARCH also has roughly 5–7% lower pinball loss than XGBoost for all five coins.

## Running it

```r
# R >= 4.3
install.packages(c("crypto2", "exuber", "rugarch", "rqPen", "xgboost", "dplyr", "tidyr",
                   "ggplot2", "patchwork", "zoo", "moments", "MASS", "Kendall", "scales"))
source("run_all.R")
```

Tables go to `output/tables/` and figures to `output/figures/`. Slow results (Monte Carlo critical values, GARCH fits and rolling forecasts, ML forecasts) are cached in `data/`. A full run from scratch takes roughly 2–3 hours on a laptop, mostly the rolling GARCH forecasts. XGBoost's quantile objective requires `xgboost` ≥ 2.0. The sample end date is fixed in `R/utils.R` for reproducibility.

## Limitations

- **Look-ahead in the bubble flag.** An episode only counts once BSADF has stayed above the critical value for log(T) days, and the critical values use the full sample length. The flag used as an ML feature therefore knows, at an episode's start, that it will last long enough to count. A real-time version would use only the BSADF statistic up to day t.
- **Few bubble days out of sample.** After the 500-day burn-in, some coins have only a handful of bubble-regime forecast days (SOL: about ten), so regime-split tests for those coins have little power.
- **Basel-style traffic light.** The regulatory rule counts exceptions of a 99% VaR over 250 days. Here it is adapted to violation rates (green ≤ α + 1pp, yellow ≤ α + 4pp) so that both levels and all sample lengths can be compared.
- **Time-series CV.** Lasso's λ is chosen by ordinary k-fold CV within each training window, not blocked CV.

## Author

Balint Tibor Nagy, MSc Quantitative Finance, WU Vienna. MIT licensed.
