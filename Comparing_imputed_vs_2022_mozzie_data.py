#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# ---------------------------------------------------------
# CLI parsing
# ---------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Linear regression comparing imputed vs 2022 mosquito egg-count differences"
    )
    parser.add_argument(
        "infile",
        help="Input CSV file (e.g. Comparing_imputed_vs_2022_mozzie_data.csv)"
    )
    return parser.parse_args()

args = parse_args()

# ---------------------------------------------------------
# Load file
# ---------------------------------------------------------

df = pd.read_csv(args.infile)

# ---------------------------------------------------------
# Sort data
# ---------------------------------------------------------

sorted_df = df.sort_values(by="ABS_DIFF_IMPUTED")

# ---------------------------------------------------------
# Convert columns to float
# ---------------------------------------------------------

x = sorted_df["ABS_DIFF_IMPUTED"].astype("float64").values
y = sorted_df["ABS_DIFF_2022_DATA"].astype("float64").values

# ---------------------------------------------------------
# Linear regression
# ---------------------------------------------------------

slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
predicted_sorted = slope * x + intercept

# ---------------------------------------------------------
# 95% confidence interval
# ---------------------------------------------------------

n = len(df)
mean_x = np.mean(x)
t_val = stats.t.ppf(0.975, df=n - 2)  # 95% CI
se = np.sqrt(np.sum((y - predicted_sorted) ** 2) / (n - 2))

conf_interval_sorted = (
    t_val
    * se
    * np.sqrt(1 / n + (x - mean_x) ** 2 / np.sum((x - mean_x) ** 2))
)

# ---------------------------------------------------------
# Plot
# ---------------------------------------------------------

plt.figure(figsize=(8, 6))

plt.scatter(
    x,
    y,
    color="black",
    label="Data points"
)

plt.plot(
    x,
    predicted_sorted,
    color="red",
    label=f"Fit line (R²={r_value**2:.2f}, p={p_value:.3f})"
)

plt.fill_between(
    x,
    predicted_sorted - conf_interval_sorted,
    predicted_sorted + conf_interval_sorted,
    color="red",
    alpha=0.2,
    label="95% CI"
)

plt.xlabel("Imputed egg count differences during In2Care intervention (2024)")
plt.ylabel("Observed egg count differences during In2Care intervention (2022)")
plt.title("Imputed vs observed egg-count differences")
plt.legend()

plt.tight_layout()
plt.show()