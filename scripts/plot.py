#!/usr/bin/env python3
"""Turn the CSV histograms written by ppmc into PNG plots.

Usage:
    python scripts/plot.py [CSV_DIR] [PNG_DIR]

Defaults: CSV_DIR = output, PNG_DIR = same as CSV_DIR. Requires matplotlib.
"""
import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")  # draw to files, no window needed
import matplotlib.pyplot as plt

# Histograms that span orders of magnitude read better on a log axis.
LOG_Y = {"dy_mass", "dy_muon_pt", "dijet_pt", "dijet_mass", "dijet_chi"}


def read_histogram(path):
    meta, lo, hi, val, err = {}, [], [], [], []
    with open(path, newline="") as f:
        rows = [line for line in f]
    for line in rows:
        if line.startswith("#"):
            text = line[1:].strip()
            if text.startswith("x:"):
                meta["x"] = text[2:].strip()
            elif text.startswith("y:"):
                meta["y"] = text[2:].strip()
            else:
                meta["title"] = text
    body = [line for line in rows if not line.startswith("#")]
    for r in csv.DictReader(body):
        lo.append(float(r["x_low"]))
        hi.append(float(r["x_high"]))
        val.append(float(r["value"]))
        err.append(float(r["error"]))
    return meta, lo, hi, val, err


def plot(path, out_dir):
    meta, lo, hi, val, err = read_histogram(path)
    centers = [(a + b) / 2 for a, b in zip(lo, hi)]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.stairs(val, lo + [hi[-1]], color="#1f5fa8", linewidth=1.6)
    ax.errorbar(centers, val, yerr=err, fmt="none", ecolor="#1f5fa8", alpha=0.6, linewidth=1)
    ax.set_title(meta.get("title", path.stem))
    ax.set_xlabel(meta.get("x", "x"))
    ax.set_ylabel(meta.get("y", "value"))
    if path.stem in LOG_Y and any(v > 0 for v in val):
        ax.set_yscale("log")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    out = out_dir / f"{path.stem}.png"
    fig.savefig(out, dpi=120)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    csv_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("output")
    out_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else csv_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    files = sorted(csv_dir.glob("*.csv"))
    if not files:
        sys.exit(f"no CSV files in {csv_dir}/ - run ppmc first")
    for path in files:
        plot(path, out_dir)


if __name__ == "__main__":
    main()
