# buruli-ulcer-mosquito-control-intervention

## Dependencies:
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
scikit-learn==1.1.2
statsmodels==0.14.6
tqdm==4.65.0
```

#### Installing dependencies with conda and pip:
```
# expected install time: less than 10 minutes

conda create -n new_env python=3.10 -c conda-forge
conda activate new_env

conda install -c conda-forge \
  epiweeks==2.3.0 \
  geopy==2.4.0 \
  joblib==1.2.0 \
  matplotlib==3.7.1 \
  numpy==1.26.4 \
  pandas==2.2.3 \
  scipy==1.11.4 \
  seaborn==0.12.2 \
  scikit-learn==1.1.2 \
  statsmodels==0.14.6 \
  tqdm==4.65.0

conda install -c conda-forge pip setuptools wheel
conda install -c conda-forge numba llvmlite
pip install --no-build-isolation pyfixest==0.40.1
```

## Make mapping plot (Fig. 1A):
```
# expected run time: less than 1 minute

# usage:
python Fig_1_v7.py
This requires the ESRI shapefile format files from the Australian Bureau of Statistics (update line 94 with path to these files):
https://www.abs.gov.au/ausstats/subscriber.nsf/log?openagent&1270055001_mb_2011_vic_shape.zip&1270.0.55.001&Data%20Cubes&85F5B2ED8E3DC957CA257801000CA953&0&July%202011&23.12.2010&Latest

```

## Spatial assignment of cases to treatment and control zones
Counts human cases occurring within 650 meters of treatment and control sites over a specified time window (defined in UNIX time), assigns each case to the nearest zone and reports treatment and control case counts and the control minus treatment counts.

```
# expected run time: less than 1 minute

# usage:
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

# expected output:
Counts:
----------------------------------------------------
Cases inside treatment zone: 1
Cases inside control zone: 6
Control-minus-treatment case count difference: 5
----------------------------------------------------
```

#### Run spatial assignment of cases to treatment and control zones across sliding windows
A wrapper Bash script is included to automate running 650m_zone_counts.py across all 136 exposure days (136 IQR windows) for each year [this is to reporduce the data in Table_S6].
```
# expected run time: less than 10 minutes

# usage:
sh run.sh 2023_cases.csv 2023_UNIX-start-end.csv 2023_report.csv
sh run.sh 2024_cases.csv 2024_UNIX-start-end.csv 2024_report.csv
sh run.sh 2025_cases.csv 2025_UNIX-start-end.csv 2025_report.csv
```

## Buruli ulcer exposure reconstruction from symptom onset dates.
Takes a case line list (symptom onset + Treatment/Control label) and a digitised incubation-period dataset (min/max ranges). It builds an empirical incubation distribution by sampling within ranges and uses Monte Carlo back-calculation (exposure = onset − incubation) to generate inferred exposure date distributions for treatment and control. Outputs include daily and weekly exposure time series, summary plots, a Gamma fit to the incubation distribution with a bootstrapped CDF band and a simple post-period Treatment/Control rate ratio with a Monte Carlo interval.
```
# expected run time: less than 5 minutes

# usage:
python buruli_exposure_backcalculation_montecarlo.py cases.csv Digitalization_IP_range.csv
```

## Peak detection for treatment and control groups
Analyses the weekly exposure time series produced by the Monte Carlo back-calculation (BU_time_series_export.csv) to estimate the timing of seasonal peaks in treatment and control arms. For each arm and year, it applies a centred rolling mean to smooth weekly counts and defines the peak as the week with the maximum smoothed value.
```
# expected run time: less than 1 minute

# usage:
python Peak_detection.py

# expected output:
Detected Peaks (Maximum Smoothed Value Only):
       Group  Year  Peak_Week  Peak_WeekNum  Peak_Count
0    Control  2023     202301             1    0.002309
1    Control  2024     202401             1    0.005219
2    Control  2025     202501             1    0.004898
3  Treatment  2023     202301             1    0.004638
4  Treatment  2024     202405             5    0.004671
5  Treatment  2025     202501             1    0.003645

Timing differences (days):
        Control  Treatment  Difference_Days_Treat_minus_Control
Year                                                           
2023 2023-01-02 2023-01-02                                    0
2024 2024-01-01 2024-01-29                                   28
2025 2024-12-30 2024-12-30                                    0
```

## Poisson rate ratio test
Takes treatment and control case counts, runs a two-sided Poisson likelihood-ratio test to compare their rates and prints the LRT statistic, p-value, incidence rate ratio (IRR) and the 95% confidence interval [this is to reporduce the data in Table_S6].

```
# expected run time: less than 1 minute

# usage:
python PRT.py [number_of_cases_in_treatment] [number_of_cases_in_control]

# example:
python PRT.py 1 6

# expected output:
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
# expected run time: less than 1 minute

# usage:
python R2_csv.py 2024_control-minus-treatment_vs_egg-counts.csv

# expected output:
R² = 0.8508
```

## Compare imputed egg counts with actual 2022 egg counts
Fits a simple linear regression comparing imputed log egg-count differences during the 2024 intervention with log egg-count differences measured during an In2Care intervention conducted in 2022.
```
# expected run time: less than 1 minute

python Comparing_imputed_vs_2022_mozzie_data.py Comparing_imputed_vs_2022_mozzie_data.csv
```


