#!/usr/bin/env python3
"""
prmtop_fix_types.py

Rename the Amber atom types of a prmtop so that CP2K reads the topology
correctly. Only the type labels (AMBER_ATOM_TYPE) change; the force-field
parameters are indexed by ATOM_TYPE_INDEX and the parameter lists, so Amber
programs give exactly the same energy with the new file.

Why this is needed (observed with CP2K 2026.2):

  1. CP2K splits the AMBER_ATOM_TYPE section by whitespace, but the section
     is written in fixed 4-character fields (20a4). A type that fills the 4
     characters (e.g. Ca2+ of the Li/Merz divalent ions) is glued to the next
     one ("Ca2+Cl-") and CP2K aborts with "Unknown element for KIND <CA2+CL->".

  2. CP2K stores the charge in the atomic kind, so it creates one kind per
     (type, charge) pair and names the extra kinds by appending a counter to
     the type name, after converting it to upper case: H -> H1, H2, H3, ...
     When such a name is also the name of another type in the topology (H1,
     H4 and H5 exist in ff14SB, C3 and C5 come from the GAFF types c3 and c5),
     atoms with different charges end up in the same kind. CP2K does not stop:
     the electrostatic energy is simply wrong (by more than 1000 kcal/mol in
     the example of the tutorial, after changing the charge of a single
     backbone H atom).

  3. Types that differ only by case (GAFF ca/ha/c/n/o and ff14SB CA/HA/C/N/O)
     become the same kind name.

The script renames (a) every type with 4 or more characters, (b) every type
whose name ends with a digit and (c) the lower-case member of every pair of
types that differ only by case. New names are built from the element symbol
followed by 'Z' and a letter (e.g. HZA, CZB), which CP2K maps to the right
element and which cannot be produced by the counter of any other type.

Usage examples:
    python3 prmtop_fix_types.py -p complex.prmtop -o AMBER.prmtop
    python3 prmtop_fix_types.py -p complex.prmtop -o AMBER.prmtop --map Ca2+:CaZ --only-long
    python3 prmtop_fix_types.py -p complex.prmtop --dry-run
"""

import argparse
import string
import sys

try:
    import parmed
    from parmed.periodic_table import Element as PT_ELEMENT
    from parmed.tools import change, outparm
except ImportError as exc:  # pragma: no cover
    sys.exit(f"ERROR: {exc}. Install ParmEd with: pip install parmed")


def parse_args():
    p = argparse.ArgumentParser(
        description="Rename Amber atom types so that CP2K reads the prmtop correctly.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-p", "--prmtop", required=True, help="Input Amber topology")
    p.add_argument("-o", "--output", default=None, help="Output Amber topology (required unless --dry-run)")
    p.add_argument("--map", nargs="+", default=[], metavar="OLD:NEW",
                   help="Explicit renaming, e.g. Ca2+:CaZ (overrides the automatic names)")
    p.add_argument("--only-long", action="store_true",
                   help="Only rename types with 4 or more characters (problem 1)")
    p.add_argument("--kinds", default=None, help="Write &KIND sections (with ELEMENT) for the new names to this file")
    p.add_argument("--map-file", default="type_renaming.dat", help="File with the OLD NEW element table")
    p.add_argument("--dry-run", action="store_true", help="Only report what would be renamed")
    a = p.parse_args()
    if not a.dry_run and not a.output:
        p.error("-o/--output is required unless --dry-run is given")
    return a


def element_of(atom):
    return PT_ELEMENT[atom.atomic_number] if atom.atomic_number > 0 else "X"


def new_name(el, taken_upper):
    """Element symbol + 'Z' + letter(s); never ends with a digit."""
    if len(el) == 1:
        cands = [el + "Z" + L for L in string.ascii_uppercase]
        cands += [el + L1 + L2 for L1 in "YXWVQ" for L2 in string.ascii_uppercase]
    else:
        cands = [el + L for L in "ZYXWVQ"]
    for cand in cands:
        if cand.upper() not in taken_upper:
            return cand
    sys.exit(f"ERROR: could not build a free type name for element {el}")


def main():
    a = parse_args()
    parm = parmed.load_file(a.prmtop)
    if not isinstance(parm, parmed.amber.AmberParm):
        sys.exit("ERROR: the topology is not an Amber prmtop")

    types, count = {}, {}
    for at in parm.atoms:
        types.setdefault(at.type, element_of(at))
        count[at.type] = count.get(at.type, 0) + 1

    reasons = {}
    for t in types:
        if len(t) >= 4:
            reasons.setdefault(t, []).append("4 characters (not readable by CP2K)")
    if not a.only_long:
        for t in types:
            if t[-1].isdigit():
                reasons.setdefault(t, []).append("ends with a digit (clashes with CP2K kind counters)")
        by_upper = {}
        for t in types:
            by_upper.setdefault(t.upper(), []).append(t)
        for up, group in by_upper.items():
            if len(group) > 1:
                # keep the upper-case (protein) type, rename the others
                keep = up if up in group else group[0]
                for t in group:
                    if t != keep:
                        reasons.setdefault(t, []).append(f"same upper-case name as {keep}")

    explicit = {}
    for item in a.map:
        if ":" not in item:
            sys.exit(f"ERROR: --map entries must be OLD:NEW, got {item}")
        old, new = item.split(":", 1)
        if len(new) > 3 or new[-1].isdigit():
            sys.exit(f"ERROR: new type name {new} must have at most 3 characters and not end with a digit")
        explicit[old] = new
        reasons.setdefault(old, []).append("explicit --map")

    taken_upper = {t.upper() for t in types if t not in reasons}
    taken_upper |= {n.upper() for n in explicit.values()}
    rename = {}
    for t in types:
        if t not in reasons:
            continue
        if t in explicit:
            rename[t] = explicit[t]
        else:
            rename[t] = new_name(types[t], taken_upper)
        taken_upper.add(rename[t].upper())

    if not rename:
        print("No atom type needs to be renamed.")
        if not a.dry_run:
            outparm(parm, a.output).execute()
            print(f"Topology copied to {a.output}")
        return

    print(f"{'old':<6s} {'new':<5s} {'element':<8s} {'atoms':>7s}  reason")
    for old, new in rename.items():
        print(f"{old:<6s} {new:<5s} {types[old]:<8s} {count[old]:>7d}  {'; '.join(reasons[old])}")
    if a.dry_run:
        return

    for old, new in rename.items():
        change(parm, "AMBER_ATOM_TYPE", f"@%{old}", new).execute()
    outparm(parm, a.output).execute()

    with open(a.map_file, "w") as fh:
        fh.write("# old new element atoms\n")
        for old, new in rename.items():
            fh.write(f"{old} {new} {types[old]} {count[old]}\n")
    if a.kinds:
        with open(a.kinds, "w") as fh:
            for old, new in rename.items():
                fh.write(f"&KIND {new.upper()}  ! was {old}\n  ELEMENT {types[old]}\n&END KIND\n")

    print(f"\n{len(rename)} type(s) renamed. Topology written to {a.output}; table in {a.map_file}")
    print("Check that an Amber single point gives the same energy with both topologies.")


if __name__ == "__main__":
    main()
