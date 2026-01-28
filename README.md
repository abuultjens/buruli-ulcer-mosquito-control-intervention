# buruli-ulcer-mosquito-control-intervention

# Dependencies:
```
epiweeks==2.3.0
geopy==2.4.0
joblib==1.2.0
matplotlib==3.7.1
numpy==1.26.4
pandas==2.2.3
pyfixest==0.40.1
scipy==1.11.4
seaborn==0.12.2
sklearn==1.1.2
statsmodels==0.14.6
tqdm==4.65.0
```

## Spatial assignment of cases to treatment and control zones
Counts human cases occurring within 650 meters of treatment and control sites over a specified time window, assigns each case to the nearest zone and reports treatment and control case counts, their difference and the total number of cases.

```
python 650m_zone_counts.py \
[cases_csv_file] \
Treatment_lat_lon.csv \
Control_lat_lon.csv \
[UNIX_start_time] \
[UNIX_end_time]

# example
python 650m_zone_counts.py \
2024_cases.csv \
Treatment_lat_lon.csv \
Control_lat_lon.csv \
1719756000 \
1725804000

Counts:
----------------------------------------------------
Cases inside treatment zone: 1
Cases inside control zone: 6
Control-minus-treatment case count difference: 5
Total unique cases in window: 31
----------------------------------------------------
```

## Buruli ulcer exposure reconstruction from symptom onset dates.
This script takes a case line list (symptom onset + Treatment/Control label) and a digitised incubation-period dataset (min/max ranges). It builds an empirical incubation distribution by sampling within ranges, then uses Monte Carlo back-calculation (exposure = onset − incubation) to generate inferred exposure date distributions for each study arm. Outputs include daily and weekly exposure time series, summary plots, a Gamma fit to the incubation distribution with a bootstrapped CDF band, and a simple post-period Treatment/Control rate ratio with a Monte Carlo interval.
```
python buruli_exposure_backcalculation_montecarlo.py cases.csv Digitalization_IP_range.csv
```

## Peak detection for treatment and control groups
This script analyses the weekly exposure time series produced by the Monte Carlo back-calculation pipeline (BU_time_series_export.csv) to estimate the timing of seasonal peaks in Treatment and Control arms. For each arm and year, it applies a centred rolling mean to smooth weekly counts and defines the peak as the week with the maximum smoothed value.
```
python Peak_detection.py
```

## Poisson rate ratio test
Takes treatment and control case counts, runs a two-sided Poisson likelihood-ratio test to compare their rates and prints the LRT statistic, p-value, incidence rate ratio (IRR) and its 95% confidence interval.

```
python PRT.py [number_of_cases_in_treatment] [number_of_cases_in_control]

# example
python PRT.py 1 6

=== Two-sided Poisson Likelihood-Ratio Test (LRT) ===
LRT statistic (D):        3.962432
Two-sided p-value:        0.046526
--------------------------------------------
λ (common under H0):      3.500000
λ_treat (H1):             1.000000
λ_control (H1):           6.000000
--------------------------------------------
IRR (treat/control):      0.166667
95% CI for IRR (score): (0.026356, 1.053927)
====================================================
```

## Regression of case-count differences on log egg-count differences
Fits a simple linear regression between the control-minus-treatment case count differences and egg count log differences and reports the coefficient of determination (R²) describing how well the linear model explains the relationship between the two variables.
```
python R2_csv.py 2024_control-minus-treatment_vs_egg-counts.csv 
R² = 0.8508
```

## Compare imputed egg counts with actual 2022 egg counts
Fits a simple linear regression comparing imputed log egg-count differences from the 2024 intervention with log egg-count differences measured during an In2Care intervention conducted in 2022.
```
python Comparing_imputed_vs_2022_mozzie_data.py Comparing_imputed_vs_2022_mozzie_data.csv
```


