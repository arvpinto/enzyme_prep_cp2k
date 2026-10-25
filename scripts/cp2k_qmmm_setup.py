#!/usr/bin/env python3
"""
cp2k_qmmm_setup.py

Prepare the &QMMM section of a CP2K input for an Amber-parameterized system
and the charge-corrected topology used with electrostatic embedding.

From an Amber topology (prmtop), the matching coordinates (rst7/inpcrd) and a
QM selection (Amber mask or a file with atom serial numbers), the script:

  1. lists the QM atoms by element as &QM_KIND sections (MM_INDEX, 1-based),
  2. finds every covalent bond cut by the QM/MM boundary and writes one
     &LINK section (IMOMM) per cut bond, with ALPHA_IMOMM estimated from the
     force-field equilibrium bond lengths, alpha = r0(QM-MM) / r0(QM-H),
  3. writes the QM &CELL from the extent of the QM atoms plus a padding,
  4. computes the nominal charge of the QM region (for &DFT CHARGE),
  5. corrects the MM charges of the residues split by the boundary so that
     the MM part of each residue carries an integer charge once the MM
     charges of the QM atoms are removed (electrostatic embedding), and
     writes the corrected topology with ParmEd (plus the equivalent ParmEd
     script, so the change can be reviewed and reproduced),
  6. writes an XYZ file of the QM region including the capping hydrogens,
  7. flags atom types that CP2K may not parse from the prmtop and writes
     the corresponding &KIND sections.

Requirements: Python >= 3.8, ParmEd, NumPy.

Example:
    python3 cp2k_qmmm_setup.py -p AMBER.prmtop -c AMBER.rst7 \\
        -m ':LIG | :44@CB,HB2,HB3,CG,ND1,HD1,CE1,HE1,NE2,CD2,HD2 | :206@CB,HB2,HB3,OG,HG' \\
        --ee-prmtop AMBER_ee.prmtop
"""

import argparse
import datetime
import math
import os
import sys
from collections import OrderedDict, defaultdict

try:
    import numpy as np
    import parmed
    from parmed.amber import AmberMask
    from parmed.periodic_table import Element as PT_ELEMENT
except ImportError as exc:  # pragma: no cover
    sys.exit(f"ERROR: {exc}. Install the requirements with: pip install parmed numpy")

# Fallback X-H equilibrium distances (angstrom) used for ALPHA_IMOMM when the
# topology has no X-H bond for the type of the QM link atom (Amber ff14SB/GAFF2 values).
FALLBACK_XH = {"C": 1.090, "N": 1.010, "O": 0.960, "S": 1.336, "P": 1.410, "Se": 1.460}
CP2K_DEFAULT_ALPHA = 1.38
VERSION = "1.0"


# --------------------------------------------------------------------------- #
# command line
# --------------------------------------------------------------------------- #
def parse_args():
    p = argparse.ArgumentParser(
        description="Generate the CP2K &QMMM include file (QM_KIND, LINK, CELL) and the "
                    "charge-corrected Amber topology for electrostatic embedding.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    io = p.add_argument_group("input")
    io.add_argument("-p", "--prmtop", required=True, help="Amber topology (prmtop)")
    io.add_argument("-c", "--coords", required=True,
                    help="Coordinates matching the topology (rst7/inpcrd/ncrst/pdb)")
    sel = io.add_mutually_exclusive_group(required=True)
    sel.add_argument("-m", "--mask", help="Amber mask selecting the QM atoms, e.g. ':LIG | :206@CB,HB2,HB3,OG,HG'")
    sel.add_argument("-f", "--index-file",
                     help="Text file with the QM atom numbers (whitespace/newline separated), e.g. "
                          "written with VMD: puts $f [$sel get serial]")
    io.add_argument("--index-base", type=int, choices=[0, 1], default=1,
                    help="Numbering of --index-file: 1 for VMD 'serial'/prmtop order, 0 for VMD 'index'")

    out = p.add_argument_group("output")
    out.add_argument("-o", "--output", default="forceeval_qmmm.cp2k.inc",
                     help="CP2K include file with &CELL, &QM_KIND and &LINK (goes inside &QMMM)")
    out.add_argument("--ee-prmtop", default=None,
                     help="Write the charge-corrected topology for electrostatic embedding to this file "
                          "(e.g. AMBER_ee.prmtop). If omitted, only the ParmEd script is written")
    out.add_argument("--parmed-script", default="charge_correction.parmed",
                     help="ParmEd script with the charge changes (run with: parmed -i <script> <prmtop>)")
    out.add_argument("--xyz", default="qm_region.xyz", help="XYZ file of the QM region with capping H atoms")
    out.add_argument("--kinds", default="kinds.cp2k.inc",
                     help="CP2K include file with &KIND sections for atom types that need an explicit element "
                          "(goes inside &SUBSYS)")
    out.add_argument("--kinds-ee", default=None,
                     help="&KIND include file for the corrected topology (default: kinds_ee.cp2k.inc, "
                          "derived from --kinds)")
    out.add_argument("--summary", default="qmmm_setup_summary.txt", help="Summary/log file")

    qm = p.add_argument_group("QM/MM options")
    qm.add_argument("--padding", type=float, default=5.0,
                    help="Distance (angstrom) added on each side of the QM atoms' extent to build the QM cell")
    qm.add_argument("--cubic", action="store_true", help="Use a cubic QM cell (largest edge in all directions)")
    qm.add_argument("--round-cell", type=float, default=0.5,
                    help="Round the QM cell edges up to a multiple of this value (angstrom, 0 = no rounding)")
    qm.add_argument("--alpha", choices=["ff", "default"], default="ff",
                    help="ff: ALPHA_IMOMM = r0(QM-MM)/r0(QM-H) from the force field; default: CP2K default (1.38)")
    qm.add_argument("--qmmm-scale", type=float, default=1.0,
                    help="QMMM_SCALE_FACTOR of the MM link atoms (scales their charge in the QM/MM potential only)")
    qm.add_argument("--charge-scheme", choices=["residue", "link", "none"], default="residue",
                    help="residue: spread the residual charge evenly over the MM atoms of each boundary residue; "
                         "link: put it on the MM link atom(s) of that residue; none: no correction")
    qm.add_argument("--fix-lj", action="store_true",
                    help="In the corrected topology, give H types with zero Lennard-Jones parameters "
                         "(e.g. HO, HW) a radius/well depth (see --lj-radius/--lj-depth), as in the CP2K biochem QM/MM howto")
    qm.add_argument("--lj-radius", type=float, default=0.3019, help="R_min/2 (angstrom) used with --fix-lj")
    qm.add_argument("--lj-depth", type=float, default=0.047, help="epsilon (kcal/mol) used with --fix-lj")
    qm.add_argument("--max-per-line", type=int, default=15,
                    help="Maximum number of indices per line in MM_INDEX (lines are continued with '\\')")
    qm.add_argument("--residual-warn", type=float, default=0.15,
                    help="Warn when the charge correction of a boundary residue exceeds this value (e)")
    return p.parse_args()


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
class Logger:
    def __init__(self, path):
        self.lines = []
        self.path = path

    def __call__(self, msg=""):
        print(msg)
        self.lines.append(msg)

    def save(self):
        with open(self.path, "w") as fh:
            fh.write("\n".join(self.lines) + "\n")


def element_symbol(atom):
    n = atom.atomic_number
    if n is None or n <= 0:
        # extra points / dummy atoms
        return "X"
    return PT_ELEMENT[n]


def residue_label(res):
    return f"{res.name}{res.idx + 1}"


def wrap_indices(indices, per_line, indent="  ", keyword="MM_INDEX"):
    """Return MM_INDEX lines using CP2K line continuation ('\\')."""
    chunks = [indices[i:i + per_line] for i in range(0, len(indices), per_line)]
    out = []
    for k, ch in enumerate(chunks):
        prefix = f"{indent}{keyword} " if k == 0 else f"{indent}{' ' * (len(keyword) + 1)}"
        cont = " \\" if k < len(chunks) - 1 else ""
        out.append(prefix + " ".join(str(i) for i in ch) + cont)
    return out


def round_up(x, step):
    if step <= 0:
        return x
    return math.ceil(x / step - 1e-9) * step


def unwrap(coords, box):
    """Unwrap coordinates around the first atom with the minimum image convention
    (orthorhombic part of the box only). Returns new coords and a flag."""
    if box is None:
        return coords, False
    L = np.array(box[:3], dtype=float)
    if np.any(L <= 0):
        return coords, False
    ref = coords[0]
    d = coords - ref
    shift = -L * np.round(d / L)
    return coords + shift, bool(np.any(np.abs(shift) > 1e-6))


def problematic_type(name):
    """Atom types CP2K cannot read from the prmtop: the AMBER_ATOM_TYPE section is
    written in fixed 4-character fields and CP2K splits it by whitespace, so a type
    with 4 characters (e.g. Ca2+) is glued to the next one."""
    return len(name) >= 4


def cp2k_kinds(parm, charges, warnings, label):
    """Reproduce the atomic kind names that CP2K builds from an Amber prmtop.

    CP2K upper-cases the atom types, groups the atoms by type name and, when a
    group has more than one charge, names its kinds NAME1, NAME2, ... in order
    of increasing charge. Returns an OrderedDict {KIND_NAME: (amber_type, element)}
    for the kinds whose name is not an element symbol, and appends warnings for
    (i) different types merged because they differ only by case and (ii) kind
    names produced by more than one type group."""
    groups = OrderedDict()
    for at in parm.atoms:
        g = groups.setdefault(at.type.upper(), {"types": OrderedDict(), "q": set(), "el": element_symbol(at)})
        g["types"].setdefault(at.type, (round(at.rmin, 6), round(at.epsilon, 6)))
        g["q"].add(float(charges[at.idx]))
    names = OrderedDict()
    owners = defaultdict(set)
    merged = []
    for up, g in groups.items():
        if len(g["types"]) > 1:
            lj = set(g["types"].values())
            merged.append(f"{'/'.join(g['types'])}{' (different LJ!)' if len(lj) > 1 else ''}")
        first_type = next(iter(g["types"]))
        if len(g["q"]) == 1:
            kind_names = [up]
        else:
            kind_names = [f"{up}{k}" for k in range(1, len(g["q"]) + 1)]
        for kn in kind_names:
            owners[kn].add(up)
            names.setdefault(kn, (first_type, g["el"]))
    if merged:
        warnings.append(f"{label}: atom types that differ only by case become one CP2K kind group: "
                        + ", ".join(merged) + ". Rename them with prmtop_fix_types.py")
    clashes = {kn: ts for kn, ts in owners.items() if len(ts) > 1}
    if clashes:
        txt = ", ".join(f"{kn} <- {'/'.join(sorted(ts))}" for kn, ts in list(clashes.items())[:8])
        more = "" if len(clashes) <= 8 else f" (+{len(clashes) - 8} more)"
        warnings.append(f"{label}: {len(clashes)} CP2K kind name(s) are produced by more than one atom type "
                        f"({txt}{more}). CP2K <= 2026.2 puts these atoms in one kind with a single charge; "
                        "rename the types with prmtop_fix_types.py")
    return OrderedDict((kn, v) for kn, v in names.items() if kn != v[1].upper())


def write_kinds(path, kind_sections, now, topname):
    with open(path, "w") as fh:
        fh.write(f"! Generated by cp2k_qmmm_setup.py v{VERSION} on {now} for {topname}\n")
        fh.write("! One &KIND section with an explicit ELEMENT for every atomic kind that CP2K builds from\n")
        fh.write("! the prmtop (one kind per atom type and charge, named TYPE or TYPE1, TYPE2, ...).\n")
        fh.write("! Include this file inside &FORCE_EVAL/&SUBSYS (use it only with that topology)\n")
        for key, (t, el) in kind_sections.items():
            fh.write(f"&KIND {key:<6s}  ! type {t}\n  ELEMENT {el}\n&END KIND\n")


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    a = parse_args()
    log = Logger(a.summary)
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")

    for f in (a.prmtop, a.coords):
        if not os.path.isfile(f):
            sys.exit(f"ERROR: file not found: {f}")

    parm = parmed.load_file(a.prmtop, xyz=a.coords)
    if not isinstance(parm, parmed.amber.AmberParm):
        sys.exit("ERROR: the topology is not an Amber prmtop")
    natom = len(parm.atoms)
    xyz = np.array(parm.coordinates, dtype=float)
    box = parm.box  # a, b, c, alpha, beta, gamma or None

    # ---------------- QM selection ----------------
    if a.mask:
        qm = sorted(AmberMask(parm, a.mask).Selected())
        sel_txt = f"Amber mask: {a.mask}"
    else:
        with open(a.index_file) as fh:
            raw = [int(tok) for tok in fh.read().split()]
        qm = sorted(set(i - a.index_base for i in raw))
        sel_txt = f"index file: {a.index_file} (base {a.index_base})"
    if not qm:
        sys.exit("ERROR: the QM selection is empty")
    if qm[0] < 0 or qm[-1] >= natom:
        sys.exit(f"ERROR: QM indices out of range (topology has {natom} atoms)")
    qmset = set(qm)

    log(f"cp2k_qmmm_setup.py v{VERSION}  ({now})")
    log(f"Topology    : {a.prmtop}  ({natom} atoms, {len(parm.residues)} residues)")
    log(f"Coordinates : {a.coords}")
    if box is not None:
        log("MM cell     : {:.4f} {:.4f} {:.4f}  angles {:.2f} {:.2f} {:.2f}".format(*box))
    log(f"QM selection: {sel_txt}")
    log(f"QM atoms    : {len(qm)}")
    log("")

    # ---------------- QM atoms per element ----------------
    kinds = OrderedDict()
    for i in qm:
        el = element_symbol(parm.atoms[i])
        if el == "X":
            sys.exit(f"ERROR: atom {i + 1} ({parm.atoms[i].name}) has no element (extra point?) and cannot be QM")
        kinds.setdefault(el, []).append(i + 1)

    # ---------------- residues in the QM region ----------------
    res_qm = OrderedDict()
    for i in qm:
        r = parm.atoms[i].residue
        res_qm.setdefault(r.idx, []).append(i)
    log("QM residues (QM atoms / residue atoms):")
    for ridx, ats in res_qm.items():
        r = parm.residues[ridx]
        names = " ".join(parm.atoms[i].name for i in ats) if len(ats) < len(r.atoms) else "(whole residue)"
        log(f"  {residue_label(r):>10s} {len(ats):4d}/{len(r.atoms):<4d} {names}")
    nwat = sum(1 for ridx in res_qm if parm.residues[ridx].name in ("WAT", "HOH", "TIP3", "SOL"))
    if nwat:
        log(f"  note: {nwat} water molecule(s) in the QM region; monitor them in MD (they can diffuse away)")
    log("")

    # ---------------- bond-type lookup ----------------
    req = {}
    for b in parm.bonds:
        if b.type is None:
            continue
        t1, t2 = b.atom1.type, b.atom2.type
        req[(t1, t2)] = req[(t2, t1)] = b.type.req

    def r0_qm_h(qatom):
        """Equilibrium QM-H distance for the type of the QM link atom."""
        # 1) an X-H bond of the same atom type exists in the topology
        cands = [v for (t1, t2), v in req.items()
                 if t1 == qatom.type and parm_types_h.get(t2, False)]
        if cands:
            return min(cands), f"{qatom.type}-H (topology)"
        el = element_symbol(qatom)
        if el in FALLBACK_XH:
            return FALLBACK_XH[el], f"{el}-H (fallback table)"
        return None, "unknown"

    parm_types_h = {}
    for at in parm.atoms:
        parm_types_h[at.type] = (at.atomic_number == 1)

    # ---------------- links ----------------
    links = []
    warnings = []
    for b in parm.bonds:
        i, j = b.atom1.idx, b.atom2.idx
        if (i in qmset) == (j in qmset):
            continue
        qi, mi = (i, j) if i in qmset else (j, i)
        qa, ma = parm.atoms[qi], parm.atoms[mi]
        r_qmmm = b.type.req if b.type is not None else float(np.linalg.norm(xyz[qi] - xyz[mi]))
        r_qh, src = r0_qm_h(qa)
        if a.alpha == "ff" and r_qh:
            alpha = r_qmmm / r_qh
        else:
            alpha = CP2K_DEFAULT_ALPHA
        links.append(dict(qm=qi, mm=mi, alpha=alpha, r_qmmm=r_qmmm, r_qh=r_qh, src=src))

        eq, em = element_symbol(qa), element_symbol(ma)
        tag = (f"{residue_label(qa.residue)}@{qa.name}(QM,{qa.type}) - "
               f"{residue_label(ma.residue)}@{ma.name}(MM,{ma.type})")
        if eq == "H" or em == "H":
            warnings.append(f"bond to hydrogen cut: {tag}. Include the H atom in the same layer as its partner")
        elif not (eq == "C" and em == "C"):
            warnings.append(f"polar bond cut ({eq}-{em}): {tag}. Prefer apolar C-C bonds (e.g. CA-CB)")
        if qa.residue is not ma.residue:
            # cut between residues: peptide bond or covalent ligand
            if {qa.name, ma.name} == {"C", "N"}:
                warnings.append(f"peptide bond cut: {tag}. Cutting at CA-C or N-CA of the neighbour keeps "
                                "the amide group intact")

    mm_count = defaultdict(int)
    qm_count = defaultdict(int)
    for L in links:
        mm_count[L["mm"]] += 1
        qm_count[L["qm"]] += 1
    for mi, n in mm_count.items():
        if n > 1:
            ma = parm.atoms[mi]
            warnings.append(f"MM atom {mi + 1} ({residue_label(ma.residue)}@{ma.name}) is bonded to {n} QM atoms; "
                            "move it into the QM region")
    for qi, n in qm_count.items():
        if n > 1:
            qa = parm.atoms[qi]
            warnings.append(f"QM atom {qi + 1} ({residue_label(qa.residue)}@{qa.name}) has {n} MM neighbours "
                            f"({n} capping H on the same atom); consider enlarging the QM region")

    log(f"QM/MM links (cut bonds): {len(links)}")
    for L in links:
        qa, ma = parm.atoms[L["qm"]], parm.atoms[L["mm"]]
        rq = f"{L['r_qh']:.3f}" if L["r_qh"] else "n/a"
        log(f"  QM {L['qm'] + 1:>7d} {residue_label(qa.residue)}@{qa.name:<4s} ({qa.type:<3s})  "
            f"MM {L['mm'] + 1:>7d} {residue_label(ma.residue)}@{ma.name:<4s} ({ma.type:<3s})  "
            f"r0(QM-MM)={L['r_qmmm']:.3f}  r0(QM-H)={rq} [{L['src']}]  ALPHA_IMOMM={L['alpha']:.3f}")
    log("")

    # ---------------- charges ----------------
    q = np.array([at.charge for at in parm.atoms], dtype=float)
    total = q.sum()
    corrections = OrderedDict()  # atom idx -> new charge
    boundary = OrderedDict()
    for ridx, ats in res_qm.items():
        r = parm.residues[ridx]
        if len(ats) == len(r.atoms):
            continue
        mm_atoms = [at.idx for at in r.atoms if at.idx not in qmset]
        boundary[ridx] = mm_atoms

    qm_charge_frac = 0.0  # nominal QM charge after correction
    log("Charge bookkeeping (e):")
    log(f"  total charge of the system           : {total:+.6f}")
    log(f"  sum of MM charges of the QM atoms    : {q[qm].sum():+.6f}")
    log("  boundary residues (MM part -> target):")
    for ridx, mm_atoms in boundary.items():
        r = parm.residues[ridx]
        q_res = sum(q[at.idx] for at in r.atoms)
        q_mm = q[mm_atoms].sum()
        target = round(q_mm)
        delta = target - q_mm
        log(f"    {residue_label(r):>10s}  residue {q_res:+.4f}  MM part {q_mm:+.4f} -> {target:+d}  "
            f"(correction {delta:+.4f} on {len(mm_atoms) if a.charge_scheme == 'residue' else 'link'} atom(s))")
        if abs(delta) > a.residual_warn:
            warnings.append(f"large charge correction ({delta:+.3f} e) in {residue_label(r)}; the boundary splits "
                            "a charged/polar group, consider moving it")
        if a.charge_scheme == "residue":
            share = delta / len(mm_atoms)
            for k in mm_atoms:
                corrections[k] = q[k] + share
        elif a.charge_scheme == "link":
            lk = [L["mm"] for L in links if parm.atoms[L["mm"]].residue.idx == ridx]
            if not lk:
                # the residue is cut, but its MM part is not the MM link atom side
                lk = mm_atoms
            share = delta / len(lk)
            for k in lk:
                corrections[k] = q[k] + share

    q_new = q.copy()
    for k, v in corrections.items():
        q_new[k] = v
    q_mm_new = sum(q_new[i] for i in range(natom) if i not in qmset)
    q_qm_nominal = total - q_mm_new
    qm_charge = int(round(q_qm_nominal))
    log(f"  MM region charge after correction    : {q_mm_new:+.6f}")
    log(f"  nominal QM charge (total - MM)       : {q_qm_nominal:+.6f}  ->  &DFT CHARGE {qm_charge}")
    if abs(q_qm_nominal - qm_charge) > 0.02:
        warnings.append(f"the QM charge ({q_qm_nominal:+.4f}) is not close to an integer; check the selection")
    if abs(total - round(total)) > 1e-3:
        warnings.append(f"the total charge of the system ({total:+.6f}) is not an integer; fix the ligand charges "
                        "before building the topology")
    log("")

    kind_sections_ee = None
    if a.ee_prmtop:
        kind_sections_ee = cp2k_kinds(parm, q_new, warnings, "corrected topology")
        if not a.kinds_ee:
            root, ext = os.path.splitext(a.kinds)
            a.kinds_ee = root.replace("kinds", "kinds_ee") if "kinds" in root else root + "_ee"
            a.kinds_ee += ext

    # ---------------- QM cell ----------------
    qxyz, split = unwrap(xyz[qm], box)
    if split:
        warnings.append("the QM region is split across the periodic boundary in these coordinates; "
                        "image/center the system (cpptraj autoimage) before running CP2K")
    # include the capping hydrogens in the extent
    caps = []
    for L in links:
        qpos = qxyz[qm.index(L["qm"])]
        v = xyz[L["mm"]] - xyz[L["qm"]]
        if box is not None:
            Lb = np.array(box[:3])
            v = v - Lb * np.round(v / Lb)
        rqh = L["r_qh"] if L["r_qh"] else 1.09
        caps.append(qpos + v / np.linalg.norm(v) * rqh)
    allpos = np.vstack([qxyz] + ([np.array(caps)] if caps else []))
    extent = allpos.max(axis=0) - allpos.min(axis=0)
    cell = extent + 2.0 * a.padding
    if a.cubic:
        cell = np.full(3, cell.max())
    cell = np.array([round_up(c, a.round_cell) for c in cell])
    log("QM cell:")
    log("  extent of QM atoms + caps            : {:.3f} {:.3f} {:.3f}".format(*extent))
    log("  QM cell (padding {:.1f} A per side)    : {:.3f} {:.3f} {:.3f}".format(a.padding, *cell))
    gmax_qm = [int(round(c)) for c in cell]
    log("  suggested QM Ewald GMAX (DFTB/SE)    : {} {} {}".format(*gmax_qm))
    if box is not None:
        gmax_mm = [int(round(c)) for c in box[:3]]
        log("  suggested MM Ewald GMAX              : {} {} {}".format(*gmax_mm))
        if np.any(cell > np.array(box[:3])):
            warnings.append("the QM cell is larger than the MM cell; reduce --padding")
    log("")

    # ---------------- CP2K atomic kinds ----------------
    # CP2K builds one atomic kind per (atom type, charge) pair. When a type has
    # several charges, the kinds are named TYPE1, TYPE2, ... in order of
    # increasing charge. The element is then guessed from the name without the
    # trailing digits, unless a &KIND section with that exact name gives it.
    # The kind names are reproduced here so that every kind gets an explicit
    # ELEMENT, and kind names shared by different atom types (a collision that
    # merges atoms with different charges in released CP2K versions up to
    # 2026.2) are reported.
    types = OrderedDict()
    for at in parm.atoms:
        types.setdefault(at.type, element_symbol(at))
    bad = OrderedDict((t, el) for t, el in types.items() if problematic_type(t))
    kind_sections = cp2k_kinds(parm, q, warnings, "original topology")

    # ---------------- write include file ----------------
    with open(a.output, "w") as fh:
        w = fh.write
        w(f"! Generated by cp2k_qmmm_setup.py v{VERSION} on {now}\n")
        w(f"! Topology: {os.path.basename(a.prmtop)}   Coordinates: {os.path.basename(a.coords)}\n")
        w(f"! QM selection: {sel_txt}\n")
        w(f"! QM atoms: {len(qm)}   QM/MM links: {len(links)}   nominal QM charge: {qm_charge}\n")
        w("! Include this file inside &FORCE_EVAL/&QMMM\n\n")
        w("! CP2K > FORCE_EVAL > QMMM > CELL\n")
        w("&CELL\n")
        w("  ABC [angstrom] {:.3f} {:.3f} {:.3f}\n".format(*cell))
        w("  ALPHA_BETA_GAMMA 90.0 90.0 90.0\n")
        w("  PERIODIC XYZ\n")
        w("&END CELL\n\n")
        w("! CP2K > FORCE_EVAL > QMMM > QM_KIND\n")
        for el, idx in kinds.items():
            w(f"&QM_KIND {el}\n")
            w("\n".join(wrap_indices(idx, a.max_per_line)) + "\n")
            w("&END QM_KIND\n")
        if links:
            w("\n! CP2K > FORCE_EVAL > QMMM > LINK\n")
        for L in links:
            qa, ma = parm.atoms[L["qm"]], parm.atoms[L["mm"]]
            w("&LINK\n")
            w(f"  ALPHA_IMOMM {L['alpha']:.3f}   ! {qa.type}(QM)-{ma.type}(MM): "
              f"{residue_label(qa.residue)}@{qa.name} - {residue_label(ma.residue)}@{ma.name}\n")
            w("  LINK_TYPE IMOMM\n")
            w(f"  MM_INDEX {L['mm'] + 1}\n")
            w(f"  QM_INDEX {L['qm'] + 1}\n")
            w(f"  QMMM_SCALE_FACTOR {a.qmmm_scale:.2f}\n")
            w("&END LINK\n")

    # ---------------- &KIND file ----------------
    write_kinds(a.kinds, kind_sections, now, os.path.basename(a.prmtop))

    # ---------------- XYZ of the QM region ----------------
    with open(a.xyz, "w") as fh:
        fh.write(f"{len(qm) + len(caps)}\n")
        fh.write(f"QM region ({len(qm)} atoms + {len(caps)} capping H), charge {qm_charge}\n")
        for k, i in enumerate(qm):
            fh.write("{:<2s} {:14.6f} {:14.6f} {:14.6f}\n".format(element_symbol(parm.atoms[i]), *qxyz[k]))
        for c in caps:
            fh.write("{:<2s} {:14.6f} {:14.6f} {:14.6f}\n".format("H", *c))

    # ---------------- charge correction ----------------
    with open(a.parmed_script, "w") as fh:
        fh.write(f"# Generated by cp2k_qmmm_setup.py v{VERSION} on {now}\n")
        fh.write(f"# Charge correction for electrostatic embedding (scheme: {a.charge_scheme})\n")
        fh.write(f"# Usage: parmed -i {os.path.basename(a.parmed_script)} {os.path.basename(a.prmtop)}\n")
        for k, v in corrections.items():
            at = parm.atoms[k]
            fh.write(f"change CHARGE @{k + 1} {v:.6f}   # {residue_label(at.residue)}@{at.name} {q[k]:+.6f} -> {v:+.6f}\n")
        if a.fix_lj:
            seen = set()
            for at in parm.atoms:
                if at.atomic_number == 1 and at.type not in seen:
                    seen.add(at.type)
                    if abs(at.epsilon) < 1e-8 or abs(at.rmin) < 1e-8:
                        fh.write(f"changeLJSingleType @%{at.type} {a.lj_radius} {a.lj_depth}\n")
        if a.ee_prmtop:
            fh.write(f"outparm {os.path.basename(a.ee_prmtop)}\n")

    if a.ee_prmtop:
        if a.charge_scheme == "none" and not a.fix_lj:
            log("Charge scheme 'none': the topology is copied without changes")
        from parmed.tools import change, changeLJSingleType, outparm
        for k, v in corrections.items():
            change(parm, "CHARGE", f"@{k + 1}", f"{v:.6f}").execute()
        if a.fix_lj:
            seen = set()
            for at in list(parm.atoms):
                if at.atomic_number == 1 and at.type not in seen:
                    seen.add(at.type)
                    if abs(at.epsilon) < 1e-8 or abs(at.rmin) < 1e-8:
                        changeLJSingleType(parm, f"@%{at.type}", str(a.lj_radius), str(a.lj_depth)).execute()
                        log(f"LJ parameters of H type {at.type} set to R_min/2={a.lj_radius} eps={a.lj_depth}")
        outparm(parm, a.ee_prmtop).execute()
        write_kinds(a.kinds_ee, kind_sections_ee, now, os.path.basename(a.ee_prmtop))
        # verify
        chk = parmed.load_file(a.ee_prmtop)
        qc = np.array([at.charge for at in chk.atoms])
        qmm = sum(qc[i] for i in range(natom) if i not in qmset)
        log(f"Corrected topology written: {a.ee_prmtop}")
        log(f"  MM region charge {qmm:+.6f}; the QM atoms keep their original MM charges "
            f"({qc[qm].sum():+.6f}), which CP2K does not use with electrostatic embedding")

    # ---------------- final report ----------------
    log("")
    log("Files written:")
    log(f"  {a.output:<32s} &CELL, &QM_KIND and &LINK sections (include inside &QMMM)")
    log(f"  {a.kinds:<32s} &KIND sections for the {len(kind_sections)} CP2K kinds of {os.path.basename(a.prmtop)} "
        "(include inside &SUBSYS)")
    if a.ee_prmtop:
        log(f"  {a.kinds_ee:<32s} &KIND sections for the {len(kind_sections_ee)} CP2K kinds of "
            f"{os.path.basename(a.ee_prmtop)}")
    log(f"  {a.xyz:<32s} QM region with capping H (visual check, standalone QM tests)")
    log(f"  {a.parmed_script:<32s} ParmEd commands for the charge correction ({len(corrections)} atoms)")
    if a.ee_prmtop:
        log(f"  {a.ee_prmtop:<32s} charge-corrected topology for ECOUPL COULOMB/GAUSS")
    log(f"  {a.summary:<32s} this summary")
    if bad:
        warnings.append("atom types with 4 characters that CP2K cannot read from the prmtop: "
                        + " ".join(bad.keys()) + ". Rename them with prmtop_fix_types.py")
    log("")
    if warnings:
        log(f"WARNINGS ({len(warnings)}):")
        for wmsg in warnings:
            log("  - " + wmsg)
    else:
        log("No warnings.")
    log("")
    log(f"Set in &FORCE_EVAL/&DFT:  CHARGE {qm_charge}   (and the MULTIPLICITY of your QM region)")
    log.save()


if __name__ == "__main__":
    main()
