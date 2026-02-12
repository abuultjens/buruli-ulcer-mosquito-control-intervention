import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from datetime import date   # for ISO week conversion

# ---------------------------------------------------------
# USER-ADJUSTABLE PARAMETERS
# ---------------------------------------------------------
SMOOTHING_WINDOW = 50
IGNORE_YEAR = 2022

# ---------------------------------------------------------
# TEXT SIZE PARAMETERS
# ---------------------------------------------------------
PEAK_LABEL_FONTSIZE = 15
AXIS_LABEL_FONTSIZE = 15
TICK_LABEL_FONTSIZE = 15
TITLE_FONTSIZE = 20
LEGEND_FONTSIZE = 15

# ---------------------------------------------------------
# LEGEND POSITION PARAMETERS
# ---------------------------------------------------------
LEGEND_X_POSITION = 0.7
LEGEND_Y_POSITION = 0.0

# ---------------------------------------------------------
# Load Data
# ---------------------------------------------------------
df = pd.read_csv("BU_time_series_export.csv")

df["Week"] = df["Week"].astype(int)
df["Year"] = df["Week"] // 100
df["WeekNum"] = df["Week"] % 100
df = df.sort_values(["Nearest zone", "Year", "Week"])

# ---------------------------------------------------------
# Smoothing Function
# ---------------------------------------------------------
def smooth_series_group_year(df, y_col="count", window=SMOOTHING_WINDOW):
    df = df.copy()
    df["smoothed"] = (
        df.groupby(["Nearest zone", "Year"])[y_col]
          .rolling(window, center=True, min_periods=1)
          .mean()
          .reset_index(level=[0,1], drop=True)
    )
    return df

df = smooth_series_group_year(df)

# ---------------------------------------------------------
# NEW Peak Detection: Maximum Smoothed Value Only
# ---------------------------------------------------------
def detect_peaks_max_only(df, ignore_year=None):
    results = []

    for (group, year), sub in df.groupby(["Nearest zone", "Year"]):
        if year == ignore_year:
            continue

        sub = sub.sort_values("Week")

        # Find index of maximum smoothed value
        idx = sub["smoothed"].idxmax()
        row = sub.loc[idx]

        results.append({
            "Group": group,
            "Year": year,
            "Peak_Week": int(row["Week"]),
            "Peak_WeekNum": int(row["WeekNum"]),
            "Peak_Count": float(row["smoothed"])
        })

    return pd.DataFrame(results)

# Run max-based peak detection
peak_results = detect_peaks_max_only(df, ignore_year=IGNORE_YEAR)

print("\nDetected Peaks (Maximum Smoothed Value Only):")
print(peak_results)

# ---------------------------------------------------------
# Convert ISO week → actual date
# ---------------------------------------------------------
def iso_week_to_date(year, week):
    return date.fromisocalendar(year, week, 1)  # Monday of ISO week

peak_results["Peak_Date"] = peak_results.apply(
    lambda r: iso_week_to_date(r["Year"], r["Peak_WeekNum"]),
    axis=1
)

# ---------------------------------------------------------
# Summaries by year, including DAY difference
# ---------------------------------------------------------
summary = (
    peak_results.pivot(index="Year", columns="Group", values="Peak_Date")
    .rename_axis(None, axis=1)
)

# Convert to pandas datetime for subtraction
for col in ["Treatment", "Control"]:
    summary[col] = pd.to_datetime(summary[col])

summary["Difference_Days_Treat_minus_Control"] = \
    (summary["Treatment"] - summary["Control"]).dt.days

print("\nTiming differences (days):")
print(summary)

# Export results
summary.to_csv("BU_peak_difference_by_year.csv", index=True)
peak_results.to_csv("BU_peak_results_complete.csv", index=False)
print("\nCSV Exported: BU_peak_difference_by_year.csv and BU_peak_results_complete.csv")

# ---------------------------------------------------------
# Plot Results
# ---------------------------------------------------------
plt.figure(figsize=(12,5))
sns.lineplot(data=df, x="Week", y="smoothed", hue="Nearest zone")

# Plot peaks
for _, r in peak_results.iterrows():
    label = f"{r['Group']} — Week {r['Peak_Week']}"
    plt.scatter(r["Peak_Week"], r["Peak_Count"], s=120, color="red")
    plt.text(
        r["Peak_Week"],
        r["Peak_Count"] - 0.0004,
        label,
        fontsize=PEAK_LABEL_FONTSIZE,
        ha="center",
        va="top"
    )

plt.title("Peak Timing (Maximum Smoothed Value)", fontsize=TITLE_FONTSIZE)
plt.xlabel("Week (YYYYWW)", fontsize=AXIS_LABEL_FONTSIZE)
plt.ylabel("Smoothed Count", fontsize=AXIS_LABEL_FONTSIZE)

plt.xticks(fontsize=TICK_LABEL_FONTSIZE)
plt.yticks(fontsize=TICK_LABEL_FONTSIZE)

plt.legend(
    title="Nearest zone",
    fontsize=LEGEND_FONTSIZE,
    title_fontsize=LEGEND_FONTSIZE,
    loc="lower center",
    bbox_to_anchor=(LEGEND_X_POSITION, LEGEND_Y_POSITION),
    ncol=2
)

plt.grid(True)
plt.tight_layout()
plt.show()
