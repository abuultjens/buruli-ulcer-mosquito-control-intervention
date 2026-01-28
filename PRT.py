#!/usr/bin/env python3
import argparse
import numpy as np
from scipy.stats import chi2
from statsmodels.stats.rates import confint_poisson_2indep

parser = argparse.ArgumentParser(
    description="Two-sided Poisson Likelihood-Ratio Test (LRT) with IRR and 95% CI."
)

parser.add_argument("treatment", type=int, help="Treatment group event count (integer >= 0)")
parser.add_argument("control", type=int, help="Control group event count (integer >= 0)")
parser.add_argument("--exposure-treatment", type=float, default=1.0, help="Treatment exposure (default 1.0)")
parser.add_argument("--exposure-control", type=float, default=1.0, help="Control exposure (default 1.0)")
parser.add_argument(
    "--ci-method",
    type=str,
    default="score",
    help="CI method for IRR from statsmodels (e.g. score, score-log). Default: score"
)

args = parser.parse_args()

# -----------------------
# Basic input checks
# -----------------------
k1 = int(args.treatment)
k2 = int(args.control)
t1 = float(args.exposure_treatment)
t2 = float(args.exposure_control)
ci_method = str(args.ci_method)

if k1 < 0 or k2 < 0:
    raise ValueError("Counts must be >= 0.")
if t1 <= 0 or t2 <= 0:
    raise ValueError("Exposures must be > 0.")

# -----------------------
# MLE rates
# -----------------------
lam0 = (k1 + k2) / (t1 + t2)   # common rate under H0
lam1 = k1 / t1                 # treatment rate under H1
lam2 = k2 / t2                 # control rate under H1

# -----------------------
# Log-likelihoods (up to additive constants)
# -----------------------
mu0_1 = lam0 * t1
mu0_2 = lam0 * t2
mu1_1 = lam1 * t1
mu2_2 = lam2 * t2

# safe log-likelihood pieces for k1 under H0
if k1 == 0:
    logL0_1 = -mu0_1
else:
    logL0_1 = k1 * np.log(mu0_1) - mu0_1

# safe log-likelihood pieces for k2 under H0
if k2 == 0:
    logL0_2 = -mu0_2
else:
    logL0_2 = k2 * np.log(mu0_2) - mu0_2

logL0 = logL0_1 + logL0_2

# safe log-likelihood pieces for k1 under H1
if k1 == 0:
    logL1_1 = -mu1_1
else:
    logL1_1 = k1 * np.log(mu1_1) - mu1_1

# safe log-likelihood pieces for k2 under H1
if k2 == 0:
    logL1_2 = -mu2_2
else:
    logL1_2 = k2 * np.log(mu2_2) - mu2_2

logL1 = logL1_1 + logL1_2

# -----------------------
# LRT statistic + p-value
# -----------------------
D = 2 * (logL1 - logL0)

# Clamp tiny negatives due to rounding
if D < 0 and D > -1e-12:
    D = 0.0

p_two = 1 - chi2.cdf(D, df=1)

# -----------------------
# IRR point estimate
# -----------------------
rate1 = k1 / t1
rate2 = k2 / t2

if rate2 == 0:
    if rate1 > 0:
        irr = np.inf
    else:
        irr = np.nan
else:
    irr = rate1 / rate2

# -----------------------
# 95% CI for IRR (score-based)
# -----------------------
ci_low, ci_high = confint_poisson_2indep(
    count1=k1,
    exposure1=t1,
    count2=k2,
    exposure2=t2,
    method=ci_method,
    compare="ratio",
    alpha=0.05
)

# -----------------------
# printing helpers
# -----------------------
def _fmt(x, digits):
    try:
        if x is None:
            return "NA"
        if isinstance(x, float) and np.isnan(x):
            return "NA"
    except Exception:
        pass
    if x == np.inf:
        return "inf"
    if x == -np.inf:
        return "-inf"
    return f"{x:.{digits}f}"

print("")
print("=== Two-sided Poisson Likelihood-Ratio Test (LRT) ===")
print("LRT statistic (D):        " + _fmt(D, 6))
print("Two-sided p-value:        " + _fmt(p_two, 6))
print("--------------------------------------------")
print("λ (common under H0):      " + _fmt(lam0, 6))
print("λ_treat (H1):             " + _fmt(lam1, 6))
print("λ_control (H1):           " + _fmt(lam2, 6))
print("--------------------------------------------")
print("IRR (treat/control):      " + _fmt(irr, 6))
print("95% CI for IRR (" + str(ci_method) + "): (" + _fmt(ci_low, 6) + ", " + _fmt(ci_high, 6) + ")")
print("====================================================")
print("")
