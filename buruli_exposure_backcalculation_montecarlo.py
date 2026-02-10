#!/usr/bin/env python3
"""
Buruli ulcer intervention analysis (script version)

This script reads a case line-list with symptom onset dates and a Treatment/Control label and a digitised incubation period
dataset (as min/max ranges) and then:
1) builds an empirical incubation period distribution by sampling within ranges
2) uses Monte Carlo sampling to infer likely exposure dates (onset - incubation)
3) aggregates inferred exposures by day and by (relative) week
4) fits a Gamma distribution to the sampled incubation periods and bootstraps a
   confidence band for the Gamma CDF
5) plots a few summaries and exports a weekly time series to CSV

Usage:
  ./buruli_intervention_analysis.py cases.csv Digitalization_IP_range.csv

Where:
  cases.csv must contain at least:
    - 'symptom onset' (UNIX seconds)
    - 'Nearest zone' (Treatment/Control)
  Digitalization_IP_range.csv must contain columns:
    - 'min'
    - 'max'
"""

import warnings
warnings.filterwarnings("ignore")

import logging
logging.getLogger().setLevel(logging.ERROR)
logging.getLogger("matplotlib").setLevel(logging.ERROR)
logging.getLogger("PIL").setLevel(logging.ERROR)

# =============================================================================
# Imports
# =============================================================================

from datetime import datetime, timezone, timedelta  # timezone/timedelta kept from notebook
import math
import itertools
import pickle
import random
import argparse
import sys
from pathlib import Path
import builtins

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import gamma, gaussian_kde
import statsmodels.api as sm
import statsmodels.formula.api as smf
import pyfixest as pf

from joblib import Parallel, delayed
from tqdm import tqdm
from epiweeks import Week

# Matplotlib log-level can be set in newer versions; ignore if not available
try:
    plt.set_loglevel("error")
except Exception:
    pass


# =============================================================================
# CLI parsing
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Buruli ulcer intervention analysis (Monte Carlo exposure back-calculation)."
    )
    parser.add_argument(
        "cases_csv",
        help="Path to case line-list CSV (must include 'symptom onset' and 'Nearest zone')."
    )
    parser.add_argument(
        "incubation_csv",
        help="Path to incubation period CSV (must include columns 'min' and 'max')."
    )
    return parser.parse_args()


args = parse_args()


# =============================================================================
# Silence ALL prints except the final RR line
# =============================================================================

_original_print = builtins.print
builtins.print = lambda *a, **k: None


# =============================================================================
# Randomness / reproducibility
# =============================================================================

# This RNG is used for the incubation-range sampling (the "empirical" distribution).
rng = np.random.default_rng(42)


# =============================================================================
# Load inputs
# =============================================================================

cases_path = Path(args.cases_csv)
incub_path = Path(args.incubation_csv)

if not cases_path.exists():
    raise FileNotFoundError(f"Cases file not found: {cases_path}")

if not incub_path.exists():
    raise FileNotFoundError(f"Incubation file not found: {incub_path}")

# Case line list (must contain 'symptom onset' and 'Nearest zone').
table_S2 = pd.read_csv(cases_path)

required_case_cols = {"symptom onset", "Nearest zone"}
missing_case = required_case_cols - set(table_S2.columns)
if missing_case:
    raise ValueError(
        f"Cases file is missing required columns: {sorted(missing_case)}. "
        f"Found columns: {list(table_S2.columns)}"
    )

# Symptom onset is stored as UNIX seconds in your file.
table_S2["Date_of_symptom_onset"] = pd.to_datetime(
    table_S2["symptom onset"], unit="s", errors="coerce"
)
if table_S2["Date_of_symptom_onset"].isna().any():
    n_bad = int(table_S2["Date_of_symptom_onset"].isna().sum())
    raise ValueError(
        f"Found {n_bad} rows where 'symptom onset' could not be parsed as UNIX seconds."
    )

# Digitised incubation period ranges (min/max, in days).
Digitalization_IP = pd.read_csv(incub_path, usecols=["min", "max"])

if not {"min", "max"}.issubset(Digitalization_IP.columns):
    raise ValueError(
        f"Incubation file must contain columns 'min' and 'max'. Found: {list(Digitalization_IP.columns)}"
    )


# =============================================================================
# Incubation-period sampling helpers
# =============================================================================

def sample_from_exact_values(df):
    """
    If you ever switch to a file of exact incubation periods (one value per row),
    this reads them from the column 'incubation_days' and returns a 1D array.
    """
    if "incubation_days" not in df.columns:
        raise ValueError("DataFrame must contain a column named 'incubation_days'.")
    return df["incubation_days"].dropna().to_numpy(dtype=float)


def sample_uniform_points_from_intervals(df, n_per_row=100, rng=None):
    """
    Turn min/max ranges into an empirical sample distribution.

    For each row i with [min_i, max_i], sample n_per_row values uniformly in
    that interval. All samples are concatenated and returned as a 1D array.
    """
    if rng is None:
        rng = np.random.default_rng()

    if not {"min", "max"}.issubset(df.columns):
        raise ValueError("DataFrame must have 'min' and 'max' columns.")

    mins = df["min"].to_numpy()
    maxs = df["max"].to_numpy()

    if np.any(maxs < mins):
        raise ValueError("Found rows where max < min.")

    widths = (maxs - mins)
    U = rng.random((len(df), n_per_row))
    samples = mins[:, None] + widths[:, None] * U
    return samples.ravel(order="C")


def build_ecdf(samples):
    """Return (x_sorted, F) for an empirical CDF."""
    x_sorted = np.sort(np.asarray(samples))
    n = x_sorted.size
    F = np.arange(1, n + 1) / n
    return x_sorted, F


def sample_from_empirical_quantile(samples, n=1, u=None, rng=None, method="linear"):
    """
    Sample from an empirical distribution using inverse-CDF sampling.

    If u is None, we draw n uniform(0,1) values and map them through np.quantile.
    """
    if rng is None:
        rng = np.random.default_rng()

    x = np.asarray(samples)
    if x.ndim != 1 or x.size == 0:
        raise ValueError("`samples` must be a non-empty 1D array.")

    if u is None:
        q = rng.random(n)
    else:
        q = np.asarray(u, dtype=float)
        if q.ndim == 0:
            q = np.repeat(q, n)
        else:
            if n != 1 and q.size != n:
                raise ValueError(f"Length of `u` ({q.size}) must equal n ({n}).")
        if np.any((q <= 0.0) | (q >= 1.0)):
            raise ValueError("All quantiles in `u` must lie strictly within (0,1).")

    try:
        draws = np.quantile(x, q, method=method)
    except TypeError:
        draws = np.quantile(x, q, interpolation=method)

    return np.atleast_1d(draws).astype(float).reshape(-1)


def plot_ecdf(samples, title="Empirical CDF of incubation period"):
    """Quick ECDF plot."""
    x_sorted, F = build_ecdf(samples)
    plt.figure()
    plt.step(x_sorted, F, where="post")
    plt.xlabel("Incubation period (days)")
    plt.ylabel("Empirical CDF")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()


def plot_pdf_hist(samples, bins="auto", title="Histogram density of incubation period"):
    """Quick histogram density plot (rough PDF)."""
    plt.figure()
    plt.hist(samples, bins=bins, density=True, edgecolor="black")
    plt.xlabel("Incubation period (days)")
    plt.ylabel("Density")
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()


def plot_ecdf_on_ax(ax, samples, label="ECDF"):
    """
    Matplotlib-3.7-compatible ECDF plotting helper.
    Replaces plt.ecdf / ax.ecdf which are only available in newer matplotlib versions.
    """
    x_sorted, F = build_ecdf(samples)
    ax.step(x_sorted, F, where="post", label=label)


# =============================================================================
# Build empirical incubation period distribution
# =============================================================================

samples = sample_uniform_points_from_intervals(Digitalization_IP, n_per_row=100, rng=rng)

# All descriptive prints are suppressed, but plots still show
plot_ecdf(samples)
plot_pdf_hist(samples, bins="auto")
plt.show()


# =============================================================================
# Monte Carlo inferred exposures (ECDF sampling)
# =============================================================================

trial_start = pd.to_datetime("2024-01-25")
trial_end = pd.to_datetime("2024-03-21")

prior_a = 1
prior_b = 1

n_iter = 200

store_incubation = []
all_counts = []
store_mu_treat = []
store_mu_control = []
store_irr = []

rng = np.random.default_rng()

for i in range(n_iter):
    draws_days = sample_from_empirical_quantile(samples, n=len(table_S2), u=None, rng=rng)
    draws_days = np.rint(draws_days).astype(int)
    store_incubation.append(draws_days)

    deltas = pd.to_timedelta(draws_days, unit="D")
    df_shift = table_S2.copy()
    df_shift["Date_of_exposure"] = df_shift["Date_of_symptom_onset"] - deltas

    counts_i = (
        df_shift.groupby(["Date_of_exposure", "Nearest zone"])
        .size()
        .reset_index(name="count")
    )
    counts_i["iter"] = i
    all_counts.append(counts_i)

    if False:
        in_window = df_shift["Date_of_exposure"].between(trial_start, trial_end, inclusive="both")
        df_win = df_shift.loc[in_window]

        k_treat = (df_win["Nearest zone"] == "Treatment").sum()
        k_control = (df_win["Nearest zone"] == "Control").sum()

        a_post_t = prior_a + k_treat
        b_post_t = prior_b + 1
        mu_treat_post_mean = a_post_t / b_post_t

        a_post_c = prior_a + k_control
        b_post_c = prior_b + 1
        mu_control_post_mean = a_post_c / b_post_c

        store_mu_treat.append(mu_treat_post_mean)
        store_mu_control.append(mu_control_post_mean)
        store_irr.append(mu_treat_post_mean / mu_control_post_mean)

df_counts = pd.concat(all_counts, ignore_index=True)


# =============================================================================
# Week-of-exposure index (relative to trial_start)
# =============================================================================

d = (df_counts["Date_of_exposure"] - trial_start).dt.days

df_counts["Week_of_exposure"] = np.where(
    d >= 0,
    (d // 7) + 1,
    -(((-d - 1) // 7) + 1)
).astype(int)

df_counts_week = (
    df_counts[["count", "iter", "Week_of_exposure", "Nearest zone"]]
    .groupby(["iter", "Week_of_exposure", "Nearest zone"], as_index=False)
    .agg(count_by_week=("count", "sum"))
)


# =============================================================================
# Average weekly exposure curves across iterations
# =============================================================================

df_counts_week_average = (
    df_counts_week
    .groupby(["Week_of_exposure", "Nearest zone"], as_index=False)
    .agg(count_by_week_mean=("count_by_week", "mean"))
)

pivot_avg = df_counts_week_average.pivot_table(
    index="Week_of_exposure",
    columns="Nearest zone",
    values="count_by_week_mean",
    aggfunc="mean"
).fillna(0)

full_weeks = np.arange(pivot_avg.index.min(), pivot_avg.index.max() + 1)
pivot_avg = pivot_avg.reindex(full_weeks, fill_value=0)
pivot_avg = pivot_avg[pivot_avg.index != 0]

x = pivot_avg.index
y_treatment = pivot_avg["Treatment"] if "Treatment" in pivot_avg else pd.Series(0, index=x)
y_control = pivot_avg["Control"] if "Control" in pivot_avg else pd.Series(0, index=x)

plt.figure(figsize=(10, 5))
plt.plot(x, y_treatment, marker="o", label="Treatment", color="C2")
plt.plot(x, y_control, marker="o", label="Control", color="C4")
plt.axvline(x=1, color="gray", linestyle="--", linewidth=1.5, label="Trial Start")
plt.axvline(x=9, color="gray", linestyle="--", linewidth=1.5, label="Trial End")
plt.title("Average weekly exposures across Monte Carlo iterations")
plt.xlabel("Week_of_exposure")
plt.ylabel("Mean exposures per week")
plt.grid(True, linestyle=":", linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()


# =============================================================================
# Sliding-window difference (Treatment - Control)
# =============================================================================

def sliding_avg_diff(pivot_avg: pd.DataFrame, n: int) -> pd.Series:
    """n-week sliding average of Treatment - Control, aligned to the first week."""
    idx = pivot_avg.index
    y_t = pivot_avg["Treatment"] if "Treatment" in pivot_avg else pd.Series(0, index=idx, dtype=float)
    y_c = pivot_avg["Control"] if "Control" in pivot_avg else pd.Series(0, index=idx, dtype=float)
    diff = (y_t - y_c).astype(float)

    kernel = np.ones(n, dtype=float) / n
    vals = np.convolve(diff.values, kernel, mode="valid")
    index = diff.index[:len(vals)]
    return pd.Series(vals, index=index, name=f"avg_diff_{n}w")


n = 8
series_n = sliding_avg_diff(pivot_avg, n=n)
series_n = series_n[series_n.index != 0]

plt.figure(figsize=(10, 5))
plt.plot(series_n.index, series_n.values, marker="o", label=f"{n}-week mean(Treatment - Control)", color="C3")
plt.axvline(x=1, color="gray", linestyle="--", linewidth=1.5, label="Trial Start")
plt.axvline(x=9, color="gray", linestyle="--", linewidth=1.5, label="Trial End")
plt.title(f"Sliding average difference ({n}-week window)")
plt.xlabel("Week_of_exposure (first week of window)")
plt.ylabel("Average difference")
plt.grid(True, linestyle=":", linewidth=0.5)
plt.legend()
plt.tight_layout()
plt.show()


# =============================================================================
# Mean differences in pre/trial/post windows (prints suppressed)
# =============================================================================

diff = (pivot_avg.get("Treatment", 0) - pivot_avg.get("Control", 0)).astype(float)
if 0 in diff.index:
    diff = diff.drop(index=0)

intervals = {
    "pre_trial_-40_to_-1": (-40, -1),
    "trial_1_to_8": (1, 8),
    "post_trial_9_to_36": (9, 36),
}

means = {}
for name, (lo, hi) in intervals.items():
    mask = (diff.index >= lo) & (diff.index <= hi)
    means[name] = diff.loc[mask].mean() if mask.any() else np.nan


# =============================================================================
# Fit a Gamma distribution to the Monte Carlo incubation draws
# =============================================================================

all_ip_draws = [x for y in store_incubation for x in y]
shape, loc, scale = gamma.fit(all_ip_draws, floc=0)

x = np.linspace(0, gamma.ppf(0.999, a=shape, scale=scale), 1000)
pdf = gamma.pdf(x, a=shape, loc=loc, scale=scale)
cdf = gamma.cdf(x, a=shape, loc=loc, scale=scale)

fig, ax = plt.subplots()
plt.hist(samples, bins=range(0, 300, 10), density=True, edgecolor="black")
plt.xlabel("Incubation period (days)")
plt.ylabel("Density")
plt.plot(x, pdf, linewidth=2, label="Gamma PDF")
plt.grid(True)
ax.text(
    x=0.5, y=1.1, s="Estimated incubation period distribution",
    fontsize=16, weight="bold", ha="center", va="bottom", transform=ax.transAxes
)
ax.text(
    x=0.5, y=1.05, s="n = 43 intervals | 200 iterations | binwidth = 10 days",
    ha="center", va="bottom", transform=ax.transAxes
)
plt.legend()
plt.show()


# =============================================================================
# Bootstrap Gamma parameters to get a CDF band (tqdm disabled cleanly)
# =============================================================================

def gamma_dist_bootstrap(data, force_loc=True):
    """Bootstrap Gamma fit by resampling and refitting."""
    random_sample = np.random.choice(data, size=len(data), replace=True)
    if force_loc:
        boot_shape, boot_loc, boot_scale = gamma.fit(random_sample, floc=0)
    else:
        boot_shape, boot_loc, boot_scale = gamma.fit(random_sample)
    return boot_shape, boot_loc, boot_scale


results = Parallel(n_jobs=10, verbose=0)(
    delayed(gamma_dist_bootstrap)(samples) for _ in tqdm(range(10000), disable=True)
)

shapes = [r[0] for r in results]
locs = [r[1] for r in results]
scales = [r[2] for r in results]

shape_0975 = np.percentile(shapes, 97.5)
shape_0025 = np.percentile(shapes, 2.5)
scale_0975 = np.percentile(scales, 97.5)
scale_0025 = np.percentile(scales, 2.5)
loc_0975 = np.percentile(locs, 97.5)
loc_0025 = np.percentile(locs, 2.5)

cdf_0975 = gamma.cdf(x, a=shape_0975, loc=loc_0975, scale=scale_0975)
cdf_0025 = gamma.cdf(x, a=shape_0025, loc=loc_0025, scale=scale_0025)

plt.figure(figsize=(8, 5))
plt.plot(x, cdf, linewidth=2, label="Gamma CDF")
plt.fill_between(x, cdf_0025, cdf_0975, color="grey", alpha=0.3)

# ---- REPLACEMENT FOR plt.ecdf(samples, ...) ----
ax = plt.gca()
plot_ecdf_on_ax(ax, samples, label="ECDF")

plt.title("Empirical CDF vs Gamma CDF (bootstrap band)")
plt.xlabel("Incubation period (days)")
plt.ylabel("Cumulative probability")
plt.grid(True)
plt.legend()
plt.show()


# =============================================================================
# Inferred exposure time series (Gamma sampling: 200 draws per case)
# =============================================================================

control_exposure = []
treat_exposure = []

for row in range(table_S2.shape[0]):
    onset = table_S2.loc[row, "Date_of_symptom_onset"]
    incubation_periods = gamma.rvs(shape, scale=scale, size=200)
    exposure_times = onset - pd.to_timedelta(incubation_periods, unit="d")

    if table_S2.loc[row, "Nearest zone"] == "Control":
        control_exposure.extend(exposure_times)
    else:
        treat_exposure.extend(exposure_times)

treat_exposure = pd.Series(treat_exposure).dt.date
control_exposure = pd.Series(control_exposure).dt.date

treat_counts = treat_exposure.value_counts().sort_index().reset_index()
treat_counts.columns = ["Date", "count"]
treat_counts["Nearest zone"] = "Treatment"

control_counts = control_exposure.value_counts().sort_index().reset_index()
control_counts.columns = ["Date", "count"]
control_counts["Nearest zone"] = "Control"

df_counts_daily = pd.concat([treat_counts, control_counts], ignore_index=True)
df_counts_daily["Date"] = pd.to_datetime(df_counts_daily["Date"])

min_date = pd.concat([pd.Series(treat_exposure), pd.Series(control_exposure)]).min()
max_date = pd.concat([pd.Series(treat_exposure), pd.Series(control_exposure)]).max()

total_daily_count = pd.MultiIndex.from_product(
    [pd.date_range(start=min_date, end=max_date, freq="D"), ["Treatment", "Control"]],
    names=["Date", "Nearest zone"]
).to_frame(index=False)

total_daily_count = (
    total_daily_count
    .merge(df_counts_daily, on=["Date", "Nearest zone"], how="left")
    .fillna({"count": 0})
)

total_daily_count["count"] = total_daily_count["count"] / (200 * table_S2.shape[0])

sns.lineplot(data=total_daily_count, x="Date", y="count", hue="Nearest zone")
plt.xlabel("Date")
plt.ylabel("Density")
plt.axvline(x=trial_start, color="gray", linestyle="--", linewidth=1.5, label="Trial Start")
plt.axvline(x=trial_end, color="gray", linestyle="--", linewidth=1.5, label="Trial End")
plt.show()

total_weekly_count = total_daily_count.copy()
total_weekly_count["Week"] = total_weekly_count["Date"].dt.strftime("%G%V")
test = total_weekly_count.groupby(["Week", "Nearest zone"])["count"].agg("sum").reset_index()

sns.lineplot(data=test, x="Week", y="count", hue="Nearest zone")
plt.xlabel("Time")
plt.ylabel("Density")
plt.xticks(
    ticks=["202101", "202201", "202301", "202401", "202501"],
    labels=["2021", "2022", "2023", "2024", "2025"]
)
plt.axvline(x="202404", color="gray", linestyle="--", linewidth=1.5, label="Trial Start")
plt.axvline(x="202412", color="gray", linestyle="--", linewidth=1.5, label="Trial End")
plt.show()

test.to_csv("BU_time_series_export.csv", index=False)


# =============================================================================
# One summary figure saved to PDF (prints suppressed)
# =============================================================================

fig = plt.figure(figsize=(12, 8))
fig.subplots_adjust(hspace=0.3, wspace=0.5)

gs = fig.add_gridspec(2, 3)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])
ax3 = fig.add_subplot(gs[:, 1:])

ax1.hist(samples, bins=range(0, 300, 10), density=True, edgecolor="black")
ax1.set_xlabel("Incubation period (days)")
ax1.set_ylabel("Density")
ax1.plot(x, pdf, linewidth=2, label="Gamma PDF")
ax1.grid(True)
ax1.text(
    0.5, 1.05, "Estimated incubation period distribution",
    fontsize=10, weight="bold", ha="center", va="bottom", transform=ax1.transAxes
)
ax1.text(
    0.5, 1.00, "n = 43 intervals | 200 iterations | binwidth = 10 days",
    fontsize=8, ha="center", va="bottom", transform=ax1.transAxes
)
ax1.legend(fontsize=8)

ax2.plot(x, cdf, linewidth=2, label="Gamma CDF")
ax2.fill_between(x, cdf_0025, cdf_0975, color="grey", alpha=0.3)

# ---- REPLACEMENT FOR ax2.ecdf(samples, ...) ----
plot_ecdf_on_ax(ax2, samples, label="ECDF")

ax2.text(
    0.5, 1.05, "Empirical CDF vs Gamma CDF",
    fontsize=10, weight="bold", ha="center", va="bottom", transform=ax2.transAxes
)
ax2.text(
    0.5, 1.00, "(bootstrap band)",
    fontsize=8, ha="center", va="bottom", transform=ax2.transAxes
)
ax2.set_xlabel("Incubation period (days)")
ax2.set_ylabel("Cumulative probability")
ax2.grid(True)
ax2.legend(fontsize=8, loc="lower right")

sns.lineplot(data=test, x="Week", y="count", hue="Nearest zone", ax=ax3)
ax3.set_xlabel("Time")
ax3.set_ylabel("Density")
ax3.text(
    0.5, 1.03, "Inferred exposure distribution (Gamma)",
    fontsize=10, weight="bold", ha="center", va="bottom", transform=ax3.transAxes
)
ax3.set_xticks(["202101", "202201", "202301", "202401", "202501"], ["2021", "2022", "2023", "2024", "2025"])
ax3.axvspan("202404", "202412", color="green", alpha=0.2, label="Intervention period")
ax3.legend()

plt.savefig("basic_incubation_period_info.pdf")


# =============================================================================
# Raw rate ratio (post period) with Monte Carlo CI
# =============================================================================

df_counts_2x2table = df_counts_week.copy()
df_counts_2x2table["post"] = df_counts_2x2table["Week_of_exposure"] >= 0
df_counts_2x2table["treat"] = df_counts_2x2table["Nearest zone"] == "Treatment"
df_counts_2x2table = (
    df_counts_2x2table
    .groupby(["post", "treat", "iter"])["count_by_week"]
    .agg("sum")
    .reset_index()
)

RR_res = []
for i in range(n_iter):
    treat_inc = df_counts_2x2table.query("post == 1 and treat == 1 and iter == @i")["count_by_week"].sum()
    ctrl_inc = df_counts_2x2table.query("post == 1 and treat == 0 and iter == @i")["count_by_week"].sum()
    rr = (treat_inc / 1) / (ctrl_inc / 1)
    RR_res.append(rr)

# =============================================================================
# FINAL AND ONLY OUTPUT
# =============================================================================

builtins.print = _original_print
print(
    f"Raw rate ratio (post) = {np.mean(RR_res):.3f} "
    f"(95% interval {np.percentile(RR_res, 2.5):.3f}–{np.percentile(RR_res, 97.5):.3f})"
)
