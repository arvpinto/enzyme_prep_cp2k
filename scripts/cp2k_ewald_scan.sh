#!/usr/bin/env bash
#
# cp2k_ewald_scan.sh - scan the Ewald ALPHA parameter of the MM part with
# single-point energies and report the real-space cutoff (RCUT) chosen by CP2K
# and the change of the energy.
#
# The template input must be a RUN_TYPE ENERGY input in which the ALPHA value
# of the MM &EWALD section is written as ALPHA_TAG and the project name as
# PROJECT_TAG, e.g.
#
#     &GLOBAL
#         RUN_TYPE ENERGY
#         PROJECT PROJECT_TAG
#     ...
#             &EWALD
#                 EWALD_TYPE SPME
#                 ALPHA ALPHA_TAG
#                 GMAX 78 83 81
#
# Usage:
#   ./cp2k_ewald_scan.sh -i energy_mm_scan.inp [-a "0.30 0.35 0.40 0.45"] [-x cp2k.psmp]
#                        [-n 4] [-t 1] [-r -113559.26] [-c 10.0] [-o ewald_scan.dat] [-k]
#
# Options:
#   -i  template input with ALPHA_TAG and PROJECT_TAG                [required]
#   -a  list of ALPHA values (angstrom^-1)       [default: "0.25 0.30 0.35 0.40 0.45 0.50"]
#   -x  CP2K executable                                         [default: cp2k.psmp]
#   -n  number of MPI tasks (0 = run without mpirun)            [default: 0]
#   -t  OpenMP threads per MPI task                             [default: 1]
#   -m  MPI launcher                                            [default: mpirun]
#   -r  reference energy in kcal/mol (e.g. the Amber/sander EPtot of the same
#       coordinates); if given, E(CP2K) - E(ref) is reported      [default: none]
#   -c  non-bonded cutoff RCUT_NB in angstrom, used to flag RCUT > RCUT_NB  [default: 10.0]
#   -o  output table                                            [default: ewald_scan.dat]
#   -k  keep the generated inputs and outputs (default: keep outputs only)
#   -h  show this help
#
# The table can be plotted with: python3 plot_ewald_scan.py -i ewald_scan.dat

set -euo pipefail

TEMPLATE=""
ALPHAS="0.25 0.30 0.35 0.40 0.45 0.50"
CP2K="cp2k.psmp"
NMPI=0
NOMP=1
MPIRUN="mpirun"
EREF=""
RCUT_NB="10.0"
OUT="ewald_scan.dat"
KEEP=0
HA2KCAL=627.509474

usage() { sed -n '2,38p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts "i:a:x:n:t:m:r:c:o:kh" opt; do
    case $opt in
        i) TEMPLATE=$OPTARG ;;
        a) ALPHAS=$OPTARG ;;
        x) CP2K=$OPTARG ;;
        n) NMPI=$OPTARG ;;
        t) NOMP=$OPTARG ;;
        m) MPIRUN=$OPTARG ;;
        r) EREF=$OPTARG ;;
        c) RCUT_NB=$OPTARG ;;
        o) OUT=$OPTARG ;;
        k) KEEP=1 ;;
        h) usage 0 ;;
        *) usage 1 ;;
    esac
done

[[ -z "$TEMPLATE" ]] && { echo "ERROR: -i is required"; usage 1; }
[[ -f "$TEMPLATE" ]] || { echo "ERROR: $TEMPLATE not found"; exit 1; }
grep -q "ALPHA_TAG" "$TEMPLATE" || { echo "ERROR: ALPHA_TAG not found in $TEMPLATE"; exit 1; }
grep -q "PROJECT_TAG" "$TEMPLATE" || { echo "ERROR: PROJECT_TAG not found in $TEMPLATE"; exit 1; }
command -v "$CP2K" > /dev/null || { echo "ERROR: $CP2K not found in PATH"; exit 1; }

export OMP_NUM_THREADS=$NOMP
RUN="$CP2K"
[[ "$NMPI" -gt 0 ]] && RUN="$MPIRUN -np $NMPI $CP2K"

printf "# %-8s %-12s %-24s %-18s %-18s %-10s\n" "ALPHA" "RCUT/A" "E/hartree" "E/kcal.mol-1" "dE_ref/kcal.mol-1" "RCUT>RCUT_NB" > "$OUT"

for A in $ALPHAS; do
    NAME="EWALD_A${A}"
    INP="$NAME.inp"
    LOG="$NAME.out"
    sed -e "s/ALPHA_TAG/$A/" -e "s/PROJECT_TAG/$NAME/" "$TEMPLATE" > "$INP"
    echo "[$(date '+%T')] ALPHA = $A"
    $RUN -i "$INP" -o "$LOG" > /dev/null 2>&1 || { echo "  WARNING: CP2K failed for ALPHA = $A (see $LOG)"; continue; }

    # the MM EWALD block is printed first; take the first RCUT reported
    RCUT=$(grep -m1 "EWALD| Real Space Cutoff" "$LOG" | awk '{print $NF}')
    E=$(grep "ENERGY| Total FORCE_EVAL" "$LOG" | tail -n 1 | awk '{print $NF}')
    [[ -z "$E" ]] && { echo "  WARNING: no energy found in $LOG"; continue; }
    EK=$(echo "$E * $HA2KCAL" | bc -l)
    DREF="nan"
    [[ -n "$EREF" ]] && DREF=$(echo "$EK - ($EREF)" | bc -l)
    FLAG=$(echo "$RCUT > $RCUT_NB" | bc -l | sed 's/1/yes/; s/0/no/')
    printf "  %-8s %-12.4f %-24.12f %-18.4f %-18s %-10s\n" "$A" "$RCUT" "$E" "$EK" \
        "$( [[ $DREF == nan ]] && echo nan || printf '%.4f' "$DREF")" "$FLAG" >> "$OUT"
    [[ $KEEP -eq 0 ]] && rm -f "$INP"
done

echo
cat "$OUT"
echo
echo "RCUT > RCUT_NB: the real-space Ewald sum is not complete within the non-bonded cutoff and the"
echo "energy is wrong (several hundred kcal/mol in the tutorial example)."
echo "Choose the smallest ALPHA giving RCUT <= RCUT_NB ($RCUT_NB A) for which the energy is stable,"
echo "and check that it does not change when GMAX (about 1 point per angstrom) is increased: a large"
echo "ALPHA needs a finer reciprocal-space grid."
