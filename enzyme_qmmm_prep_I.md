<br>
<h2> <p align="center"> <b>I - Setup of the enzyme-ligand complex with AmberTools</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
PDB structure  ->  pdb4amber / PROPKA  ->  antechamber / parmchk2  ->  tleap  ->  complex.prmtop + complex.rst7
</pre>

---

<br/>

<p align="justify">In this section we build the Amber topology (prmtop) and coordinates (rst7) of the solvated enzyme-ligand complex. Every decision taken here (protonation states, ligand charges, ions, box) is carried into the QM/MM model, so it pays off to make them carefully and to write them down.</p>

<br/>

<h3> <p align="center"><b> Cleaning the PDB structure </b></p> </h3>

<br/>

- <p align="justify">Start by inspecting the structure in VMD and the header of the PDB file: missing residues and atoms (REMARK 465 and 470), alternate conformations, chain breaks, crystallization additives, metal ions and cofactors. Keep the crystallographic waters, especially those in the active site, and remove buffer molecules and other additives that are not part of your model. Missing loops have to be modelled (e.g. with MODELLER) before going further.</p>

- <p align="justify">Insertion codes (residues 36A, 60B, 184A, ...) are common in proteases and antibodies, which follow the numbering of a reference protein. pdb4amber handles them poorly: on the thrombin structure it inserted 46 TER records at residues with insertion codes, which would break the chain in tleap. We therefore renumber the residues sequentially first, with <a href="http://www.bonvinlab.org/pdb-tools/" target="_blank">pdb-tools</a> (pip install pdb-tools), blank the insertion code column and keep a map between both numberings:</p>

<pre style="color: white; background-color: black;">
grep -E "^(ATOM|HETATM|TER|END)" thrombin_2ZFF.pdb | pdb_reres -1 | \
    sed -E 's/^((ATOM  |HETATM).{20})./\1 /' > thrombin_seq.pdb

paste &lt;(awk '/^ATOM|^HETATM/{print substr($0,18,10)}' thrombin_2ZFF.pdb | uniq) \
      &lt;(awk '/^ATOM|^HETATM/{print substr($0,18,10)}' thrombin_seq.pdb | uniq) > residue_map.txt

grep -E "HIS H  57 |ASP H 102 |GLY H 193 |SER H 195 " residue_map.txt
HIS H  57   HIS H  44
ASP H 102   ASP H 100
GLY H 193   GLY H 204
SER H 195   SER H 206
</pre>

<br/>

- <p align="justify">pdb4amber renames residues to their Amber names, removes the hydrogen atoms (they are added back by tleap with the correct names), detects disulfide bonds (CYS -> CYX) and writes TER records between chains.</p>

<pre style="color: white; background-color: black;">
pdb4amber -i thrombin_seq.pdb -o protein_noH.pdb --nohyd
</pre>

<pre style="color: white; background-color: black;">
---------- Missing heavy atom(s)

None

cat protein_noH_sslink
29 45
120 265
174 188
202 232
</pre>

<br/>

- <p align="justify">protein_noH_sslink lists the disulfide bonds (pairs of residue numbers), protein_noH_renum.txt maps the input numbering to the output one and protein_noH_nonprot.pdb holds everything that is not protein. Keep residue_map.txt at hand: the catalytic residues of thrombin are His57, Asp102 and Ser195 in the literature, but His44, Asp100 and Ser206 in the topology.</p>

- <p align="justify">Check the TER records (here, 2: one per chain) and every chain break reported. The version of pdb4amber used here prints many "gap of ... A between residues i and i+2" lines for this file that are not real gaps (all consecutive C-N distances are about 1.33 &Aring;), so always confirm a reported gap in VMD. If the structure has a real gap that you do not want to model, cap both ends (ACE and NME) instead of letting tleap join two distant residues with a long peptide bond. In the example, the two chains of thrombin were already capped.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Protonation states </b></p> </h3>

<br/>

- <p align="justify">Estimate the pKa values of the titratable residues with PROPKA (or H++, PDB2PQR) at the pH you want to simulate:</p>

<pre style="color: white; background-color: black;">
propka3 thrombin_2ZFF.pdb
</pre>

<pre style="color: white; background-color: black;">
SUMMARY OF THIS PREDICTION
       Group      pKa  model-pKa
   ...
   ASP 102 H     3.61       3.80
   ASP 189 H     6.92       3.80
   ...
   HIS  57 H     6.94       6.50
   ...
</pre>

<br/>

- <p align="justify">The predicted values are a starting point, not the answer. Inspect the hydrogen bond network of every residue with a pKa close to the working pH and of every residue involved in catalysis. For the catalytic triad of a serine protease, His57 donates a hydrogen bond from N&delta;1 to Asp102 and accepts one from Ser195 O&gamma; at N&epsilon;2, so it has to be neutral and protonated at N&delta;1 (HID). Asp189, at the bottom of the S1 pocket, has an elevated pKa (6.9) but stays deprotonated at pH 7. The protonation states of the residues that take part in the reaction must be consistent with the mechanism you want to study (for some mechanisms you will need to build more than one model).</p>

- <p align="justify">The protonation states are set by renaming the residues in the PDB file: HID, HIE or HIP for histidines, ASH for neutral aspartate, GLH for neutral glutamate, LYN for neutral lysine and CYX for cysteines in disulfide bonds (already done by pdb4amber). The default HIS is read by tleap as HIE.</p>

<pre style="color: white; background-color: black;">
sed -e 's/HIS H  44/HID H  44/' \
    -e 's/HIS H  67/HIE H  67/' -e 's/HIS H  88/HIE H  88/' \
    -e 's/HIS H 117/HIE H 117/' -e 's/HIS H 243/HIE H 243/' protein_noH.pdb | grep -v CONECT > protein_amber.pdb
</pre>

<br/>

- <p align="justify">Note the grep -v CONECT: pdb4amber writes CONECT records for the disulfide bonds and tleap reads them. If you also define the disulfides with the bond command, as we do below, tleap aborts when saving the topology with "1-4: cannot add bond ... This may be caused by duplicate bond specifications". Use one of the two, not both.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Ligand parameters </b></p> </h3>

<br/>

- <p align="justify">The ligand needs all its hydrogen atoms and the right protonation state before the charges are derived. The D-Phe-Pro-benzylamide inhibitor has a primary amine that is protonated at pH 7 (net charge +1). Here we use GAFF2 atom types and AM1-BCC charges:</p>

<pre style="color: white; background-color: black;">
antechamber -i lig_5.sdf -fi mdl -o LIG.mol2 -fo mol2 -c bcc -nc 1 -at gaff2 -rn LIG -pf y
</pre>

<br/>

- <p align="justify">antechamber writes the charges with 6 decimals and the rounding often leaves a small residual. Check the sum before going further, since a non-integer total charge is carried into every Ewald sum of the MM and QM/MM calculations:</p>

<pre style="color: white; background-color: black;">
awk '/@&lt;TRIPOS&gt;ATOM/{f=1;next}/@&lt;TRIPOS&gt;BOND/{f=0}f{s+=$9}END{printf "%.6f\n",s}' LIG.mol2
0.997001
</pre>

<pre style="color: white; background-color: black;">
python3 mol2_charge_fix.py -i LIG.mol2 -o LIG_fixed.mol2 -c 1
Atoms               : 52
Original charge sum : +0.997001
Target charge       : +1
Corrected sum       : +1.000000
Atoms modified      : 52
Written             : LIG_fixed.mol2
</pre>

<br/>

- <p align="justify">mol2_charge_fix.py spreads the residual over the atoms in steps of 10<sup>-6</sup> e (option --method largest puts it on the atom with the largest |q| instead). Then generate the missing force-field parameters and check the penalty scores in the frcmod file (high penalties deserve a closer look and, possibly, a reparameterization):</p>

<pre style="color: white; background-color: black;">
parmchk2 -i LIG_fixed.mol2 -f mol2 -o LIG.frcmod -s gaff2
</pre>

<br/>

- <p align="justify">For production work, RESP charges (HF/6-31G* on optimized geometries, consistent with the ff14SB charge model) are preferred to AM1-BCC. Metal ions bound to the protein need a dedicated parameterization (e.g. MCPB.py), and covalently bound ligands have to be built as non-standard residues. Remember that the QM region removes the MM description of whatever it contains, but the MM parameters of the ligand still drive the MD of Section II and the environment of every QM/MM calculation.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Building the system with tleap </b></p> </h3>

<br/>

- <p align="justify">We load the ff14SB protein force field, GAFF2 and the TIP3P water model with its Joung-Cheatham/Li-Merz ion parameters, combine protein and ligand, define the disulfide bonds, solvate the complex in an orthorhombic box with 12 &Aring; of water around the solute, neutralize it and add NaCl to about 0.15 M. CP2K reads the orthorhombic cell directly into &CELL, which keeps the QM/MM setup simple.</p>

<pre style="color: white; background-color: black;">
# tleap.in
source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

loadamberparams LIG.frcmod
LIG = loadmol2 LIG_fixed.mol2
prot = loadpdb protein_amber.pdb

bond prot.29.SG prot.45.SG
bond prot.120.SG prot.265.SG
bond prot.174.SG prot.188.SG
bond prot.202.SG prot.232.SG

complex = combine {prot LIG}
check complex
charge complex

solvatebox complex TIP3PBOX 12.0
addions complex Cl- 0
addionsrand complex Na+ NIONS Cl- NIONS
charge complex

saveamberparm complex complex.prmtop complex.rst7
savepdb complex complex.pdb
quit
</pre>

<br/>

- <p align="justify">The number of ion pairs for 0.15 M can be estimated from the number of water molecules added by solvatebox (run tleap once without the addions lines and read the "Added ... residues" line), N = 0.15 x N<sub>water</sub> / 55.5:</p>

<pre style="color: white; background-color: black;">
NWAT=$(grep "Added" tleap_dry.log | awk '{print $2}')
NIONS=$(python3 -c "print(round(0.15*$NWAT/55.5))")
sed "s/NIONS/$NIONS/g" tleap.in > tleap_run.in
tleap -f tleap_run.in > tleap.log
</pre>

<pre style="color: white; background-color: black;">
Total unperturbed charge:   5.000000
...
  Added 12638 residues.
...
Total unperturbed charge:  -0.000000
Exiting LEaP: Errors = 0; Warnings = 2; Notes = 0.
</pre>

<br/>

- <p align="justify">Read the whole tleap.log. The two warnings above are the expected non-zero charge of the unit before neutralization. Look for missing parameters, atoms added that should have been in the PDB, close contacts and unexpected bonds. Open complex.pdb in VMD and confirm the protonation states, the disulfides and the position of the ligand.</p>

- <p align="justify">tleap writes the solvent at the end of the topology, after the ions, so the crystallographic waters do not keep their position in the residue list. In the example the ligand is residue 291 and the ions follow it. From now on, take every residue and atom number from complex.pdb or from the prmtop, never from the original PDB file.</p>

<br/>

<p align="justify">The input files of this section are in <a href="AMBER_INPUTS.zip">AMBER_INPUTS.zip</a>.</p>

<br/>
