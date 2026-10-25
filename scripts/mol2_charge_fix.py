#!/usr/bin/env python3
"""
mol2_charge_fix.py

Make the partial charges of a Tripos MOL2 file (e.g. the AM1-BCC output of
antechamber) sum exactly to an integer net charge.

antechamber writes charges with 6 decimals, and the rounding often leaves a
small residual (e.g. +1.001000 instead of +1). tleap then builds a system
with a non-integer total charge, which is carried into the Ewald sums of the
MM and QM/MM calculations. This script removes the residual by spreading it
over the atoms (default) or by placing it on the atom(s) with the largest
absolute charge.

Usage examples:
    python3 mol2_charge_fix.py -i LIG.mol2 -o LIG_fixed.mol2 -c 1
    python3 mol2_charge_fix.py -i LIG.mol2 -o LIG_fixed.mol2 -c -1 --method largest
"""

import argparse
import sys
from decimal import Decimal, ROUND_HALF_UP


def parse_args():
    p = argparse.ArgumentParser(
        description="Correct the partial charges of a MOL2 file so that they sum to an integer.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-i", "--input", required=True, help="Input MOL2 file")
    p.add_argument("-o", "--output", required=True, help="Output MOL2 file")
    p.add_argument("-c", "--charge", type=int, default=None,
                   help="Target net charge (default: nearest integer to the current sum)")
    p.add_argument("-m", "--method", choices=["spread", "largest"], default="spread",
                   help="spread: distribute the residual over all atoms in steps of 1e-6; "
                        "largest: put the residual on the atom with the largest |q|")
    p.add_argument("--decimals", type=int, default=6, help="Number of decimals written for the charges")
    p.add_argument("-q", "--quiet", action="store_true", help="Do not print the summary")
    return p.parse_args()


def main():
    a = parse_args()
    with open(a.input) as fh:
        lines = fh.readlines()

    # locate the ATOM block
    try:
        start = next(i for i, l in enumerate(lines) if l.strip().upper() == "@<TRIPOS>ATOM") + 1
    except StopIteration:
        sys.exit("ERROR: no @<TRIPOS>ATOM section found in " + a.input)
    end = next((i for i in range(start, len(lines)) if lines[i].startswith("@<TRIPOS>")), len(lines))

    atom_idx = [i for i in range(start, end) if lines[i].strip()]
    quantum = Decimal(1).scaleb(-a.decimals)  # e.g. 0.000001
    charges = []
    for i in atom_idx:
        f = lines[i].split()
        if len(f) < 9:
            sys.exit("ERROR: atom line without a charge column:\n" + lines[i])
        charges.append(Decimal(f[8]).quantize(quantum, rounding=ROUND_HALF_UP))

    total = sum(charges)
    target = Decimal(a.charge) if a.charge is not None else total.quantize(Decimal(1), rounding=ROUND_HALF_UP)
    residual = total - target
    nsteps = int((residual / quantum).to_integral_value(rounding=ROUND_HALF_UP))

    new = list(charges)
    if nsteps != 0:
        step = -quantum if nsteps > 0 else quantum
        if a.method == "largest":
            order = sorted(range(len(new)), key=lambda k: abs(new[k]), reverse=True)
            new[order[0]] += step * abs(nsteps)
        else:
            # spread one quantum at a time, starting from the atoms with the largest |q|
            order = sorted(range(len(new)), key=lambda k: abs(new[k]), reverse=True)
            for n in range(abs(nsteps)):
                new[order[n % len(order)]] += step

    # rewrite the charge column, keeping the rest of each line untouched
    for k, i in enumerate(atom_idx):
        f = lines[i].split()
        old_q = f[8]
        line = lines[i].rstrip("\n")
        pos = line.rfind(old_q)
        newq = f"{new[k]:.{a.decimals}f}"
        width = max(len(old_q), len(newq))
        line = line[:pos].rstrip() + " " + newq.rjust(width) + line[pos + len(old_q):]
        lines[i] = line + "\n"

    with open(a.output, "w") as fh:
        fh.writelines(lines)

    if not a.quiet:
        print(f"Atoms               : {len(new)}")
        print(f"Original charge sum : {total:+.{a.decimals}f}")
        print(f"Target charge       : {int(target):+d}")
        print(f"Corrected sum       : {sum(new):+.{a.decimals}f}")
        print(f"Atoms modified      : {len([1 for x, y in zip(charges, new) if x != y])}")
        print(f"Written             : {a.output}")


if __name__ == "__main__":
    main()
