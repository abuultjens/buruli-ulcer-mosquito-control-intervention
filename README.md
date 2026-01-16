# buruli-ulcer-mosquito-control-intervention





```
python 650m_zone_counts.py \
../66-MEDIAN_Inner_northwest_2024_cases_symptom_UPDATE-LAT-LON.csv \
../Treatment_lat_lon.csv \
../Control_lat_lon.csv \
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


```
python PRT_USED_v2.py 1 6

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
