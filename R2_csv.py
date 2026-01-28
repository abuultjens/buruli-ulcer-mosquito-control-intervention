#!/usr/bin/env python3
import pandas as pd
import sys
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

def calc_r2(csv_file: str) -> float:
    # Read comma-separated 2-column CSV file
    df = pd.read_csv(csv_file, sep=",", header=None, usecols=[0, 1])
    df = df.dropna()

    if df.shape[0] == 0:
        raise ValueError("No valid data found in the file.")

    x = df.iloc[:, 0].values.reshape(-1, 1)  # predictor
    y = df.iloc[:, 1].values                 # response

    # Fit linear regression
    model = LinearRegression()
    model.fit(x, y)
    y_pred = model.predict(x)

    # Calculate R²
    r2 = r2_score(y, y_pred)
    return r2

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python calc_r2.py <file.csv>")
        sys.exit(1)

    csv_file = sys.argv[1]
    try:
        r2 = calc_r2(csv_file)
        print(f"R² = {r2:.4f}")
    except Exception as e:
        print(f"Error: {e}")
