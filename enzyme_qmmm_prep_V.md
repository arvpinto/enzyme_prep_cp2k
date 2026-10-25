<br>
<h2> <p align="center"> <b>V - Charge correction for electrostatic embedding</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
AMBER.prmtop  ->  charge_correction.parmed (ParmEd)  ->  AMBER_ee.prmtop  ->  &QMMM ECOUPL COULOMB / GAUSS
</pre>

---

<br/>

<h3> <p align="center"><b> Why the MM charges have to be corrected </b></p> </h3>

<br/>

- <p align="justify">With mechanical embedding (ECOUPL NONE), the QM and MM regions interact through the force field: the QM atoms keep their MM point charges and Lennard-Jones parameters for the QM/MM terms. With electrostatic embedding (ECOUPL COULOMB for DFTB and semiempirical methods, ECOUPL GAUSS for DFT), the QM electron density feels the MM point charges directly, and the MM charges of the QM atoms are removed: the QM atoms are represented by their nuclei and their density.</p>

- <p align="justify">Amber residues have integer charges, but their fragments do not. When the boundary cuts the C&alpha;-C&beta; bond, the side chain goes to the QM region and the backbone atoms left in the MM region carry a fractional charge. For the three cut residues of the example, the MM charges of the 74 QM atoms add up to +0.2044 e, so removing them leaves the MM region with a net charge of -0.2044 e, and puts a spurious fractional charge right next to the QM region:</p>

<pre style="color: white; background-color: black;">
Charge bookkeeping (e):
  total charge of the system           : -0.000000
  sum of MM charges of the QM atoms    : +0.204400
  boundary residues (MM part -> target):
         HID44  residue -0.0000  MM part -0.0075 -> +0  (correction +0.0075 on 6 atom(s))
        ASP100  residue -1.0000  MM part -0.1419 -> +0  (correction +0.1419 on 6 atom(s))
        SER206  residue -0.0000  MM part -0.0550 -> +0  (correction +0.0550 on 6 atom(s))
  MM region charge after correction    : -0.000000
  nominal QM charge (total - MM)       : -0.000000  ->  &DFT CHARGE 0
</pre>

<br/>

- <p align="justify">The correction restores an integer charge in the MM part of every residue split by the boundary. Following the <a href="https://manual.cp2k.org/trunk/methods/qm_mm/builtin.html" target="_blank">CP2K biochemical QM/MM example</a>, the residual of each residue is spread evenly over its MM atoms, here the six backbone atoms N, H, CA, HA, C and O. The integer charge that the QM part then carries (ligand +1, Asp102 side chain -1, the other two side chains 0) gives the charge of the QM region, &DFT CHARGE 0.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> The corrected topology </b></p> </h3>

<br/>

- <p align="justify">cp2k_qmmm_setup.py writes the changes as a ParmEd script, and applies it when --ee-prmtop is given. Each line keeps the residue, atom name and old and new charges as a comment, so that the correction can be reviewed and reproduced:</p>

<pre style="color: white; background-color: black;">
cat charge_correction.parmed
# Charge correction for electrostatic embedding (scheme: residue)
# Usage: parmed -i charge_correction.parmed AMBER.prmtop
change CHARGE @650 -0.414450   # HID44@N -0.415700 -> -0.414450
change CHARGE @651 0.273150   # HID44@H +0.271900 -> +0.273150
change CHARGE @652 0.020050   # HID44@CA +0.018800 -> +0.020050
change CHARGE @653 0.089350   # HID44@HA +0.088100 -> +0.089350
change CHARGE @665 0.598550   # HID44@C +0.597300 -> +0.598550
change CHARGE @666 -0.566650   # HID44@O -0.567900 -> -0.566650
change CHARGE @1641 -0.492650   # ASP100@N -0.516300 -> -0.492650
...
change CHARGE @3280 -0.558733   # SER206@O -0.567900 -> -0.558733
outparm AMBER_ee.prmtop
</pre>

<pre style="color: white; background-color: black;">
parmed -i charge_correction.parmed AMBER.prmtop     # same result as --ee-prmtop
</pre>

<br/>

- <p align="justify">The corrected topology is used only for electrostatic embedding, in both PARM_FILE_NAME and CONN_FILE_NAME. Calculations with mechanical embedding and pure MM calculations keep the original AMBER.prmtop. The QM atoms keep their original charges in AMBER_ee.prmtop, so the total charge of this file is not an integer (+0.2044 in the example); that is expected, since CP2K does not use them.</p>

<pre style="color: white; background-color: black;">
&FORCE_EVAL
    METHOD QMMM
    ...
    &MM
        &FORCEFIELD
            ...
            PARM_FILE_NAME AMBER_ee.prmtop
        ...
    &QMMM
        ECOUPL COULOMB
        @INCLUDE forceeval_qmmm.cp2k.inc
    &END QMMM
    &SUBSYS
        &TOPOLOGY
            CONN_FILE_FORMAT AMBER
            CONN_FILE_NAME AMBER_ee.prmtop
            ...
        &END TOPOLOGY
        @INCLUDE kinds_ee.cp2k.inc
    &END SUBSYS
</pre>

<br/>

- <p align="justify">Note the kinds_ee.cp2k.inc file. CP2K builds one atomic kind per atom type and charge, so changing charges changes the list of kinds, and each topology needs its own &KIND file (Section IV).</p>

- <p align="justify">The output of CP2K shows the charge of the classical system after the QM atoms are removed, which is the quickest check that the correction was applied. With electrostatic embedding and the corrected topology:</p>

<pre style="color: white; background-color: black;">
 FORCEFIELD| Total charge of the classical system:                  0.0000018932
</pre>

<p align="justify">and with the uncorrected topology:</p>

<pre style="color: white; background-color: black;">
 FORCEFIELD| Total charge of the classical system:                 -0.2044001141
</pre>

<br/>

---

<br/>

<h3> <p align="center"><b> Options </b></p> </h3>

<br/>

- <p align="justify">--charge-scheme link puts the whole residual of each residue on its MM link atom (C&alpha;) instead of spreading it. It keeps the rest of the residue unchanged, but concentrates the correction on the atom closest to the QM region, which can overpolarize it. --charge-scheme none writes no correction.</p>

- <p align="justify">The script warns when the correction of a residue is larger than 0.15 e (--residual-warn): it means the boundary splits a polar or charged group, and a different cut is a better solution than a large correction.</p>

- <p align="justify">--fix-lj gives a Lennard-Jones radius and well depth (R<sub>min</sub>/2 = 0.3019 &Aring;, &epsilon; = 0.047 kcal/mol, the values of a GAFF2 hydroxyl hydrogen used in the CP2K example) to the hydrogen types that have none in the force field, such as HO (hydroxyl) and HW (TIP3P water). Without repulsion, these positive MM charges can get too close to the QM density in electrostatic embedding QM/MM MD. The change applies to every atom of those types, including all water molecules, so it modifies the MM model: use it when you observe the problem, and consistently in all calculations.</p>

- <p align="justify">QMMM_SCALE_FACTOR (Section IV) is complementary: it scales the charge of the MM link atom in the QM/MM potential only, and does not change the classical charge balance discussed here.</p>

<br/>
