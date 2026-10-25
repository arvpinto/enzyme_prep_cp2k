#!/usr/bin/env python3
"""
plot_ewald_scan.py

Plot the table written by cp2k_ewald_scan.sh: the real-space cutoff chosen by
CP2K (RCUT) and the MM energy as a function of the Ewald ALPHA parameter.
The figure is shown on screen (unless --no-show) and saved to a file.

Usage:
    python3 plot_ewald_scan.py -i ewald_scan.dat -o ewald_scan.png --rcut-nb 10
"""

import argparse
import sys

import numpy as np

try:
    import matplotlib
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover
    sys.exit("ERROR: matplotlib is required (pip install matplotlib)")


def parse_args():
    p = argparse.ArgumentParser(description="Plot the Ewald ALPHA scan of cp2k_ewald_scan.sh",
                                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--input", default="ewald_scan.dat", help="Table written by cp2k_ewald_scan.sh")
    p.add_argument("-o", "--output", default="ewald_scan.png", help="Output image")
    p.add_argument("--dpi", type=int, default=600, help="Resolution of the saved image")
    p.add_argument("--rcut-nb", type=float, default=10.0, help="Non-bonded cutoff RCUT_NB (angstrom), drawn as a line")
    p.add_argument("--relative", action="store_true",
                   help="Plot the energy relative to the reference column (E_CP2K - E_ref) when available")
    p.add_argument("--no-show", action="store_true", help="Do not open a window, only save the image")
    return p.parse_args()


def main():
    a = parse_args()
    if a.no_show:
        matplotlib.use("Agg")
    rows = []
    with open(a.input) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            rows.append((float(f[0]), float(f[1]), float(f[3]), f[4]))
    if not rows:
        sys.exit(f"ERROR: no data in {a.input}")
    rows.sort()
    alpha = np.array([r[0] for r in rows])
    rcut = np.array([r[1] for r in rows])
    ek = np.array([r[2] for r in rows])
    has_ref = all(r[3] != "nan" for r in rows)
    if a.relative and has_ref:
        energy = np.array([float(r[3]) for r in rows])
        elabel = r"$E_\mathrm{CP2K} - E_\mathrm{ref}$ / kcal mol$^{-1}$"
    else:
        energy = ek - ek[-1]
        elabel = r"$E - E(\alpha_\mathrm{max})$ / kcal mol$^{-1}$"

    fig, ax1 = plt.subplots(figsize=(6.0, 4.0))
    c1, c2 = "#6676FC", "#FF4545"
    ax1.plot(alpha, rcut, "o-", color=c1, lw=1.8, ms=6, label="RCUT")
    ax1.axhline(a.rcut_nb, color=c1, ls="--", lw=1.0)
    ax1.text(alpha[-1], a.rcut_nb, f"RCUT_NB = {a.rcut_nb:g} " + r"$\mathrm{\AA}$", color=c1,
             va="bottom", ha="right", fontsize=9)
    ax1.set_xlabel(r"Ewald ALPHA / $\mathrm{\AA}^{-1}$")
    ax1.set_ylabel(r"RCUT / $\mathrm{\AA}$", color=c1)
    ax1.tick_params(axis="y", colors=c1)

    ax2 = ax1.twinx()
    ax2.plot(alpha, energy, "s-", color=c2, lw=1.8, ms=6, label="Energy")
    ax2.set_ylabel(elabel, color=c2)
    ax2.tick_params(axis="y", colors=c2)

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="lower center", frameon=False, fontsize=9)
    for ax in (ax1, ax2):
        for s in ("top",):
            ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(a.output, dpi=a.dpi)
    print(f"Figure saved to {a.output} ({a.dpi} dpi)")
    if not a.no_show:
        plt.show()


if __name__ == "__main__":
    main()
