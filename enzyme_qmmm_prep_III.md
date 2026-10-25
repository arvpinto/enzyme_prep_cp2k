<br>
<h2> <p align="center"> <b>III - Selection of the QM region</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
AMBER.prmtop + AMBER.rst7  ->  QM selection (Amber mask or VMD serials)  ->  qm_atoms.dat
</pre>

---

<br/>

<p align="justify">The QM region is the most consequential choice of a QM/MM study: it decides which chemistry the model can describe, how the environment polarizes the reacting atoms, and how much each calculation costs. There is no universal recipe, but there are well established rules, and the setup script of the next section checks the ones that can be checked automatically.</p>

<br/>

<h3> <p align="center"><b> Guidelines </b></p> </h3>

<br/>

- <p align="justify">Include every atom that forms or breaks a bond, changes its protonation state or its oxidation state along the reaction, together with the groups conjugated with them (aromatic rings, carboxylates, amides).</p>

- <p align="justify">Include the groups that interact strongly and specifically with the reacting atoms and change along the reaction: hydrogen bond donors that stabilize a developing negative charge (e.g. the oxyanion hole), metal ions and their first coordination shell, catalytic waters.</p>

- <p align="justify">Cut covalent bonds only between apolar sp<sup>3</sup> carbon atoms, preferably the C&alpha;-C&beta; bond of side chains or the C&alpha;-C bond of the backbone. Never cut polar bonds (C-N, C-O) or peptide bonds, and keep the boundary at least three bonds away from any bond that changes during the reaction. Avoid a QM atom with two MM neighbours (two capping hydrogens on the same atom) and an MM atom bonded to two QM atoms.</p>

- <p align="justify">Make sure the QM method you will use has parameters for every element in the region (basis sets and pseudopotentials for DFT, Slater-Koster tables for DFTB, etc.).</p>

- <p align="justify">Larger is not always better within a fixed budget: a larger QM region with a cheaper method, or less sampling, can be worse than a smaller, well chosen one. The properties of interest (barriers, reaction energies, charges) should be tested for convergence with respect to the size of the QM region, or at least for the effect of adding the next shell of residues.</p>

- <p align="justify">Water molecules can be part of the QM region, but in QM/MM MD they diffuse away. For geometry optimizations, NEB and dimer calculations (as in the TS tutorial) this is not a problem.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> QM region of the example </b></p> </h3>

<br/>

<p align="justify">For thrombin we take the whole ligand and the side chains of the catalytic triad, His57 (44), Asp102 (100) and Ser195 (206), cutting the three C&alpha;-C&beta; bonds. The Amber mask selection is compact: the whole ligand, plus everything in the three residues that is not a backbone atom.</p>

<pre style="color: white; background-color: black;">
:LIG | (:44,100,206 &amp; !@N,H,CA,HA,C,O)
</pre>

<br/>

- <p align="justify">The same selection in VMD, loading the topology and the snapshot (VMD resid follows the sequential numbering of the prmtop), written to a file with the atom serial numbers (1-based, as CP2K and the prmtop):</p>

<pre style="color: white; background-color: black;">
vmd -parm7 AMBER.prmtop -rst7 AMBER.rst7
</pre>

<pre style="color: white; background-color: black;">
set sel [atomselect top "resname LIG or (resid 44 100 206 and not name N H CA HA C O)"]
$sel num
74
set outfile [open qm_atoms.dat w]
puts $outfile [$sel get serial]
close $outfile
</pre>

<br/>

- <p align="justify">VMD is the most convenient way to build selections based on distances, e.g. all residues with an atom within 4 &Aring; of the ligand (same residue as within 4 of resname LIG). Selections like this one cut every residue at its peptide bonds, so they have to be edited by hand into side chains or whole backbone units before they make a sensible QM region.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Including backbone groups </b></p> </h3>

<br/>

<p align="justify">In serine proteases the backbone N-H groups of Gly193 and Ser195 form the oxyanion hole, which stabilizes the tetrahedral intermediates of acylation and deacylation. If they have to be in the QM region, each amide group must be included whole, cutting at the neighbouring C&alpha;-C bonds. Including only the N-H atoms of Gly193 (204) gives two polar cuts and a QM nitrogen with two capping hydrogens, which the setup script of the next section reports:</p>

<pre style="color: white; background-color: black;">
python3 cp2k_qmmm_setup.py -p AMBER.prmtop -c AMBER.rst7 -m ':LIG | (:44,100,206 &amp; !@N,H,CA,HA,C,O) | :204@N,H'
...
WARNINGS (4):
  - polar bond cut (N-C): GLY204@N(QM,N) - GLY204@CA(MM,CX). Prefer apolar C-C bonds (e.g. CA-CB)
  - polar bond cut (N-C): GLY204@N(QM,N) - GLU203@C(MM,C). Prefer apolar C-C bonds (e.g. CA-CB)
  - peptide bond cut: GLY204@N(QM,N) - GLU203@C(MM,C). Cutting at CA-C or N-CA of the neighbour keeps the amide group intact
  - QM atom 3251 (GLY204@N) has 2 MM neighbours (2 capping H on the same atom); consider enlarging the QM region
</pre>

<br/>

<p align="justify">The correct selection takes the C=O of the preceding residue together with the N-H (and the C&alpha; of glycine with its hydrogens), so that every cut is a C&alpha;-C bond. This extended QM region has 87 atoms and 6 links, and passes all checks:</p>

<pre style="color: white; background-color: black;">
:LIG | (:44,100 &amp; !@N,H,CA,HA,C,O) | :203@C,O | :204@N,H,CA,HA2,HA3 | :205@C,O | :206@N,H,CA,HA,CB,HB2,HB3,OG,HG
</pre>

<pre style="color: white; background-color: black;">
QM/MM links (cut bonds): 6
  QM    3272 SER206@CA   (CX )  MM    3279 SER206@C    (C  )  ...  ALPHA_IMOMM=1.396
  QM    3268 ASP205@C    (C  )  MM    3260 ASP205@CA   (CX )  ...  ALPHA_IMOMM=1.396
  QM    3253 GLY204@CA   (CX )  MM    3256 GLY204@C    (C  )  ...  ALPHA_IMOMM=1.396
  QM    3249 GLU203@C    (C  )  MM    3238 GLU203@CA   (CX )  ...  ALPHA_IMOMM=1.396
  QM    1645 ASP100@CB   (2C )  MM    1643 ASP100@CA   (CX )  ...  ALPHA_IMOMM=1.400
  QM     654 HID44@CB   (CT )  MM     652 HID44@CA   (CX )  ...  ALPHA_IMOMM=1.400
...
No warnings.
</pre>

<br/>

<p align="justify">In the rest of the tutorial we use the 74-atom QM region (ligand and catalytic triad side chains), which is enough to show every step. Whatever the selection, write it down exactly (mask or VMD selection) together with the snapshot it was applied to: the atom numbers in the CP2K input are only meaningful for that topology.</p>

<br/>
