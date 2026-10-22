"""
Project Fulcrum - fetch the input data.

Downloads antibody_agonism.csv from the public copy in this repository and
checks it is the panel the analysis expects. Run this before 01_analysis.py.
Exits non-zero if the file is missing or does not match.
"""
import os, sys, urllib.request

URL = "https://raw.githubusercontent.com/meganghuang/dataset/main/project_fulcrum/antibody_agonism.csv"
DEST = os.path.join(os.path.dirname(os.path.abspath(__file__)), "antibody_agonism.csv")

EXPECTED_COLS = ["ab_id", "site", "epitope_dist_nm", "log10_kon", "log10_koff",
                 "log10_KD", "receptor_density", "plate", "agonism_pct"]
EXPECTED_ROWS = 500  # column order in the CSV is irrelevant; presence is what matters

if not os.path.exists(DEST):
    print("downloading", URL)
    urllib.request.urlretrieve(URL, DEST)
else:
    print("already present:", DEST)

import pandas as pd
df = pd.read_csv(DEST)
problems = []
missing = [c for c in EXPECTED_COLS if c not in df.columns]
if missing:
    problems.append(f"missing columns {missing}; got {list(df.columns)}")
if len(df) != EXPECTED_ROWS:
    problems.append(f"{len(df)} rows, expected {EXPECTED_ROWS}")
if df.isna().sum().sum():
    problems.append("missing values present; the panel should be complete")
# log10_KD is defined as log10_koff - log10_kon; the analysis relies on it
dev = (df.log10_KD - (df.log10_koff - df.log10_kon)).abs().max()
if dev > 1e-9:
    problems.append(f"log10_KD != log10_koff - log10_kon (max deviation {dev:.2g})")

if problems:
    print("INPUT CHECK FAILED:")
    for p in problems:
        print("  -", p)
    sys.exit(1)

print(f"ok: {len(df)} rows, {(df.agonism_pct == 0).sum()} exact zeros in the outcome, "
      f"log10_KD identity holds to {dev:.1e}")
