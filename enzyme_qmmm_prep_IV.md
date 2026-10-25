<br>
<h2> <p align="center"> <b>IV - Building the &QMMM section of the CP2K input</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
&FORCE_EVAL
    METHOD QMMM
    ...
    &QMMM
        ECOUPL ...
        @INCLUDE forceeval_qmmm.cp2k.inc
    &END QMMM
...
</pre>

---

<br/>

<p align="justify">The &QMMM section tells CP2K which atoms are QM (&QM_KIND), how the covalent bonds cut by the boundary are capped (&LINK), and the size of the box in which the QM problem is solved (&CELL). With large QM regions these sections are long lists of atom numbers, and errors in them are easy to make and hard to spot. The cp2k_qmmm_setup.py script writes them from the topology, the coordinates and the QM selection, together with the other files that the QM/MM inputs need.</p>

<br/>

<h3> <p align="center"><b> Atom types that CP2K can read </b></p> </h3>

<br/>

<p align="justify">Before any CP2K calculation, the atom types of the topology have to be checked. CP2K (versions up to 2026.2, the latest release when this tutorial was written) has three problems with Amber atom type names, and only the first one stops the run:</p>

- <p align="justify">CP2K splits the AMBER_ATOM_TYPE section of the prmtop by whitespace, but the section is written in fixed 4-character fields. A type that fills the 4 characters, such as Ca2+ of the divalent ion parameters, is glued to the next one and the run aborts with "Unknown element for KIND &lt;CA2+CL-&gt;".</p>

- <p align="justify">CP2K keeps the charge in the atomic kind, so it builds one kind per atom type and charge, and names the extra kinds by appending a counter to the upper-case type name: the backbone amide hydrogens (type H, with a few different charges) become H1, H2, H3, ... When one of these names is also the name of another type (ff14SB has H1, H4 and H5; N gives N3, which is a type of lysine), the atoms of both types are put in the same kind, with a single charge. CP2K does not stop or warn. In the example, the original topology happens to give almost the right energy, but changing the charge of one backbone H atom (which is exactly what the charge correction of Section V does) moved the CP2K MM energy by 1148 kcal/mol, while Amber changed by less than 1 kcal/mol. This was fixed in the development version of CP2K in October 2026 (<a href="https://github.com/cp2k/cp2k/pull/6191" target="_blank">cp2k/cp2k#6191</a>).</p>

- <p align="justify">Types that differ only by case (GAFF2 ca, ha, c, n, o and ff14SB CA, HA, C, N, O) are merged into the same kind name.</p>

<p align="justify">prmtop_fix_types.py renames the offending types: every type with 4 or more characters, every type whose name ends with a digit, and the lower-case member of every pair that differs only by case. The new names are the element followed by Z and a letter, which no counter can produce. Only the labels change: the parameters are indexed by number in the prmtop.</p>

<pre style="color: white; background-color: black;">
python3 prmtop_fix_types.py -p complex.prmtop -o AMBER.prmtop
</pre>

<pre style="color: white; background-color: black;">
old    new   element    atoms  reason
H1     HZA   H            479  ends with a digit (clashes with CP2K kind counters)
O2     OZA   O             76  ends with a digit (clashes with CP2K kind counters)
H4     HZB   H             14  ends with a digit (clashes with CP2K kind counters)
C8     CZA   C            148  ends with a digit (clashes with CP2K kind counters)
N2     NZA   N             60  ends with a digit (clashes with CP2K kind counters)
N3     NZB   N             22  ends with a digit (clashes with CP2K kind counters)
H5     HZC   H              5  ends with a digit (clashes with CP2K kind counters)
ca     CZB   C             12  same upper-case name as CA
c3     CZC   C              3  ends with a digit (clashes with CP2K kind counters)
c      CZD   C              2  same upper-case name as C
o      OZB   O              2  same upper-case name as O
n      NZC   N              1  same upper-case name as N
c5     CZE   C              4  ends with a digit (clashes with CP2K kind counters)
ha     HZD   H             10  same upper-case name as HA
hc     HZE   H              6  same upper-case name as HC
h1     HZF   H              5  ends with a digit (clashes with CP2K kind counters); same upper-case name as H1

16 type(s) renamed. Topology written to AMBER.prmtop; table in type_renaming.dat
</pre>

<br/>

<p align="justify">Confirm that Amber gives exactly the same energy with both files (sp_sander.in, in <a href="CP2K_INPUTS.zip">CP2K_INPUTS.zip</a>, is a single point without the long-range dispersion correction, which CP2K does not apply):</p>

<pre style="color: white; background-color: black;">
sander -O -i sp_sander.in -o sp_complex.out -p complex.prmtop -c AMBER.rst7 -r tmp.rst7
sander -O -i sp_sander.in -o sp_AMBER.out   -p AMBER.prmtop   -c AMBER.rst7 -r tmp.rst7
grep -A5 "NSTEP       ENERGY" sp_complex.out sp_AMBER.out
</pre>

<pre style="color: white; background-color: black;">
 BOND    =      748.6821  ANGLE   =     1982.1971  DIHED      =     3569.3896
 VDWAALS =    20286.1551  EEL     =  -179221.8245  HBOND      =        0.0000
 1-4 VDW =     1028.9822  1-4 EEL =    11129.5421  RESTRAINT  =        0.0000
</pre>

<br/>

---

<br/>

<h3> <p align="center"><b> Running the setup script </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
python3 cp2k_qmmm_setup.py -p AMBER.prmtop -c AMBER.rst7 \
    -m ':LIG | (:44,100,206 &amp; !@N,H,CA,HA,C,O)' --ee-prmtop AMBER_ee.prmtop
</pre>

<p align="justify">or, with the file of atom serial numbers written in VMD:</p>

<pre style="color: white; background-color: black;">
python3 cp2k_qmmm_setup.py -p AMBER.prmtop -c AMBER.rst7 -f qm_atoms.dat --ee-prmtop AMBER_ee.prmtop
</pre>

<pre style="color: white; background-color: black;">
cp2k_qmmm_setup.py v1.0
Topology    : AMBER.prmtop  (42455 atoms, 12932 residues)
Coordinates : AMBER.rst7
MM cell     : 77.6793 82.5830 81.2552  angles 90.00 90.00 90.00
QM selection: Amber mask: :LIG | (:44,100,206 &amp; !@N,H,CA,HA,C,O)
QM atoms    : 74

QM residues (QM atoms / residue atoms):
       HID44   11/17   CB HB2 HB3 CG ND1 HD1 CE1 HE1 NE2 CD2 HD2
      ASP100    6/12   CB HB2 HB3 CG OD1 OD2
      SER206    5/11   CB HB2 HB3 OG HG
      LIG291   52/52   (whole residue)

QM/MM links (cut bonds): 3
  QM    3274 SER206@CB   (2C )  MM    3272 SER206@CA   (CX )  r0(QM-MM)=1.526  r0(QM-H)=1.090 [2C-H (topology)]  ALPHA_IMOMM=1.400
  QM    1645 ASP100@CB   (2C )  MM    1643 ASP100@CA   (CX )  r0(QM-MM)=1.526  r0(QM-H)=1.090 [2C-H (topology)]  ALPHA_IMOMM=1.400
  QM     654 HID44@CB   (CT )  MM     652 HID44@CA   (CX )  r0(QM-MM)=1.526  r0(QM-H)=1.090 [CT-H (topology)]  ALPHA_IMOMM=1.400

Charge bookkeeping (e):
  ...
  nominal QM charge (total - MM)       : -0.000000  ->  &DFT CHARGE 0

QM cell:
  extent of QM atoms + caps            : 12.371 12.570 10.562
  QM cell (padding 5.0 A per side)    : 22.500 23.000 21.000
  suggested QM Ewald GMAX (DFTB/SE)    : 22 23 21
  suggested MM Ewald GMAX              : 78 83 81
...
Files written:
  forceeval_qmmm.cp2k.inc          &CELL, &QM_KIND and &LINK sections (include inside &QMMM)
  kinds.cp2k.inc                   &KIND sections for the 246 CP2K kinds of AMBER.prmtop (include inside &SUBSYS)
  kinds_ee.cp2k.inc                &KIND sections for the 262 CP2K kinds of AMBER_ee.prmtop
  qm_region.xyz                    QM region with capping H (visual check, standalone QM tests)
  charge_correction.parmed         ParmEd commands for the charge correction (18 atoms)
  AMBER_ee.prmtop                  charge-corrected topology for ECOUPL COULOMB/GAUSS
  qmmm_setup_summary.txt           this summary

No warnings.

Set in &FORCE_EVAL/&DFT:  CHARGE 0   (and the MULTIPLICITY of your QM region)
</pre>

<br/>

<p align="justify">The charge bookkeeping and the charge-corrected topology are the subject of Section V. Here we go through the include file and the other outputs.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> QM atoms </b></p> </h3>

<br/>

- <p align="justify">One &QM_KIND section is written per element, listing the 1-based indices of the QM atoms in the topology (MM_INDEX). Long lists are split over several lines with the CP2K line continuation character (\), which avoids the problems of very long input lines:</p>

<pre style="color: white; background-color: black;">
! CP2K > FORCE_EVAL > QMMM > QM_KIND
&QM_KIND C
  MM_INDEX 654 657 660 663 1645 1648 3274 4627 4628 4629 4630 4631 4632 4633 4634 \
           4635 4638 4639 4640 4641 4642 4645 4646 4647 4648 4649 4650 4651
&END QM_KIND
&QM_KIND H
  MM_INDEX 655 656 659 661 664 1646 1647 3275 3276 3278 4653 4654 4655 4656 4657 \
           4658 4659 4660 4661 4662 4663 4664 4665 4666 4667 4668 4669 4670 4671 4672 \
           4673 4674 4675 4676 4677 4678
&END QM_KIND
&QM_KIND N
  MM_INDEX 658 662 4637 4644 4652
&END QM_KIND
&QM_KIND O
  MM_INDEX 1649 1650 3277 4636 4643
&END QM_KIND
</pre>

<br/>

- <p align="justify">The element of each QM atom comes from the atomic number stored in the prmtop, not from the atom type, so GAFF types such as ca or ns are handled correctly. When the QM method is DFT, the &SUBSYS section needs one &KIND section per QM element with its basis set and pseudopotential.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> QM/MM links </b></p> </h3>

<br/>

- <p align="justify">Every bond between a QM and an MM atom is found from the topology and gets a &LINK section. With the IMOMM scheme, CP2K caps the QM atom with a hydrogen placed along the QM-MM bond, at a distance that is the QM-MM distance divided by ALPHA_IMOMM, and projects the forces on the capping hydrogen onto the QM and MM atoms. The CP2K default (1.38) corresponds to a C-C bond; the script estimates ALPHA_IMOMM from the force-field equilibrium distances, alpha = r<sub>0</sub>(QM-MM) / r<sub>0</sub>(QM-H), using an X-H bond of the same atom type found in the topology, or a table of typical X-H distances (C 1.09, N 1.01, O 0.96, S 1.336 &Aring;) otherwise. For a C&alpha;-C&beta; cut, 1.526/1.090 = 1.400. Use --alpha default to keep 1.38 everywhere.</p>

<pre style="color: white; background-color: black;">
! CP2K > FORCE_EVAL > QMMM > LINK
&LINK
  ALPHA_IMOMM 1.400   ! 2C(QM)-CX(MM): SER206@CB - SER206@CA
  LINK_TYPE IMOMM
  MM_INDEX 3272
  QM_INDEX 3274
  QMMM_SCALE_FACTOR 1.00
&END LINK
...
</pre>

<br/>

- <p align="justify">QMMM_SCALE_FACTOR scales the charge of the MM link atom in the QM/MM electrostatic potential only (the classical MM electrostatics are not changed). The MM link atom sits about 1.5 &Aring; from a QM atom and its point charge can overpolarize the QM density; values between 0 and 1 (e.g. 0.5 to 0.8) are often used to attenuate it. The default of the script is 1.0 (no scaling), which can be changed for all links with --qmmm-scale. Whatever you choose, keep it fixed for every calculation that will be compared.</p>

- <p align="justify">The script reports a warning when the selection cuts a polar bond, a peptide bond or a bond to hydrogen, when an MM atom is bonded to more than one QM atom or a QM atom to more than one MM atom (Section III shows an example).</p>

<br/>

---

<br/>

<h3> <p align="center"><b> QM cell </b></p> </h3>

<br/>

<pre style="color: white; background-color: black;">
! CP2K > FORCE_EVAL > QMMM > CELL
&CELL
  ABC [angstrom] 22.500 23.000 21.000
  ALPHA_BETA_GAMMA 90.0 90.0 90.0
  PERIODIC XYZ
&END CELL
</pre>

<br/>

- <p align="justify">The QM cell is built from the extent of the QM atoms and the capping hydrogens plus a padding on each side (--padding, 5 &Aring; by default), rounded up to 0.5 &Aring;. CP2K centers the QM atoms in this cell (CENTER EVERY_STEP is the default), and warns if a QM atom gets too close to its border ("One or few QM atoms are within the SKIN of the quantum box"). For plane-wave based methods (GPW/GAPW) the cell has to contain the electron density, and the energy should be converged with respect to its size: increase the padding until the energy no longer changes appreciably. For DFTB and semiempirical methods the cell only defines the periodic Ewald treatment of the QM charges, and the Ewald GMAX of the &DFT/&POISSON section follows the "1 point per angstrom" rule (22 23 21 here).</p>

- <p align="justify">The script unwraps the QM atoms across the periodic boundary before measuring them, and warns if the coordinates had the QM region split by the boundary (image the system again with cpptraj, Section II).</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Atom kinds </b></p> </h3>

<br/>

- <p align="justify">CP2K builds one atomic kind per atom type and charge (TYPE, or TYPE1, TYPE2, ... in order of increasing charge, as explained at the top of this section) and needs the element of each kind. Unless a &KIND section with that exact name gives the ELEMENT, it takes the name without the trailing digits and looks it up in its built-in lists of Amber, CHARMM and GROMOS types, or in the periodic table. Most GAFF2 types and the renamed types are not in those lists, and with the GAFF2 ligand the first CP2K test (an MM single point) aborted with:</p>

<pre style="color: white; background-color: black;">
 * [ABORT]        Unknown element for KIND &lt;NS&gt;. This problem can be fixed     *
 *  \___/  specifying properly elements in PDB or specifying a KIND section or *
</pre>

<p align="justify">The setup script reproduces the kind names that CP2K builds and writes one &KIND section with an explicit ELEMENT for each of them. Since the names depend on the charges, there is one file per topology: kinds.cp2k.inc for AMBER.prmtop (246 kinds in the example) and kinds_ee.cp2k.inc for AMBER_ee.prmtop (262 kinds, because the charge correction creates new ones). Each goes inside &SUBSYS of the inputs that use the corresponding topology:</p>

<pre style="color: white; background-color: black;">
&KIND HC1     ! type HC
  ELEMENT H
&END KIND
&KIND HC2     ! type HC
  ELEMENT H
&END KIND
...
&KIND NS      ! type ns
  ELEMENT N
&END KIND
...
</pre>

<pre style="color: white; background-color: black;">
&FORCE_EVAL
    ...
    &SUBSYS
        ...
        @INCLUDE kinds.cp2k.inc
    &END SUBSYS
</pre>

<br/>

- <p align="justify">While building the kind names, the script also checks for the problems described at the top of this section (4-character types, kind names produced by more than one type, types that differ only by case) and reports them as warnings. With the topology written by prmtop_fix_types.py there are none.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Checking the QM region </b></p> </h3>

<br/>

- <p align="justify">qm_region.xyz contains the QM atoms and the capping hydrogens at their IMOMM positions (along the QM-MM bond, at the QM-H equilibrium distance). Open it in VMD over the full system and confirm that the region is what you intended, that no capping hydrogen points into another QM atom, and that the region looks chemically sensible: count the electrons, check the formal charge (&DFT CHARGE) and multiplicity.</p>

<pre style="color: white; background-color: black;">
head -n 3 qm_region.xyz
77
QM region (74 atoms + 3 capping H), charge 0
C       47.685844      32.849030      39.401936
</pre>

<pre style="color: white; background-color: black;">
vmd -parm7 AMBER.prmtop -rst7 AMBER.rst7 qm_region.xyz
</pre>

<br/>

- <p align="justify">All options of the script can be listed with python3 cp2k_qmmm_setup.py -h. The summary printed on screen is also saved to qmmm_setup_summary.txt; keep it with your inputs.</p>

<br/>
