import pandas as pd
import numpy as np
from geopy.distance import geodesic
import sys

# Set the zone radius (in km)
ZONE_RADIUS_KM = 0.65

# calculate distances (km)
def haversine_distances(sample_coords, site_coords):
    return np.array([
        [geodesic((lat1, lon1), (lat2, lon2)).km for lat2, lon2 in site_coords]
        for lat1, lon1 in sample_coords
    ])

def main(cases_file, treatment_sites_file, control_sites_file, start_unix, end_unix):
    # Load data
    cases = pd.read_csv(cases_file)
    treatment_sites = pd.read_csv(treatment_sites_file)
    control_sites = pd.read_csv(control_sites_file)

    # Convert unix_time column to numeric
    cases["unix_time"] = pd.to_numeric(cases["unix_time"], errors="coerce")

    # Filter cases within a specified time window
    cases_window = cases[
        (cases["unix_time"] >= start_unix) &
        (cases["unix_time"] <= end_unix)
    ].copy()

    # If no cases, still print zeros
    if cases_window.empty:
        print("Counts:")
        print("----------------------------------------------------")
        print("Cases Inside Treatment Zone:", 0)
        print("Cases Inside Control Zone:", 0)
        print("Control-minus-treatment case count difference:", 0)
        print("----------------------------------------------------")
        return

    # Convert lat/lon to numpy arrays
    treatment_coords = treatment_sites[["lat", "lon"]].to_numpy()
    control_coords = control_sites[["lat", "lon"]].to_numpy()

    # Compute distances to treatment and control sites
    treatment_distances = haversine_distances(
        cases_window[["lat", "lon"]].to_numpy(), treatment_coords
    )
    control_distances = haversine_distances(
        cases_window[["lat", "lon"]].to_numpy(), control_coords
    )

    # Assign each case to the nearest zone
    cases_window["nearest_zone"] = np.where(
        treatment_distances.min(axis=1) < control_distances.min(axis=1),
        "Treatment",
        "Control"
    )

    # Determine if cases are within the zone radius
    cases_window["within_treatment_zone"] = (
        treatment_distances.min(axis=1) <= ZONE_RADIUS_KM
    ).astype(int)

    cases_window["within_control_zone"] = (
        control_distances.min(axis=1) <= ZONE_RADIUS_KM
    ).astype(int)

    # Count cases inside each zone
    treatment_case_count = int(
        cases_window[cases_window["nearest_zone"] == "Treatment"]["within_treatment_zone"].sum()
    )

    control_case_count = int(
        cases_window[cases_window["nearest_zone"] == "Control"]["within_control_zone"].sum()
    )

    # Control minus treatment difference
    control_minus_treatment = control_case_count - treatment_case_count

    # Total cases in window
    total_cases_in_window = int(cases_window.shape[0])

    # Print results
    print("Counts:")
    print("----------------------------------------------------")
    print("Cases inside treatment zone:", treatment_case_count)
    print("Cases inside control zone:", control_case_count)
    print("Control-minus-treatment case count difference:", control_minus_treatment)
    print("Total unique cases in window:", total_cases_in_window)
    print("----------------------------------------------------")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(
            "Usage: python script.py <cases_file> <treatment_sites_file> "
            "<control_sites_file> <start_unix> <end_unix>"
        )
        sys.exit(1)

    cases_file = sys.argv[1]
    treatment_sites_file = sys.argv[2]
    control_sites_file = sys.argv[3]
    start_unix = int(sys.argv[4])
    end_unix = int(sys.argv[5])

    main(cases_file, treatment_sites_file, control_sites_file, start_unix, end_unix)
