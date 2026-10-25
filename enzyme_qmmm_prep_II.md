<br>
<h2> <p align="center"> <b>II - MD equilibration and choice of the starting structure</b> </p></h2>

<br/>

<pre style="color: white; background-color: black;">
min1 -> min2 -> heat (NVT) -> npt -> prod  ->  cpptraj (analysis, clustering, imaging)  ->  snapshot.rst7
</pre>

---

<br/>

<p align="justify">A QM/MM study should start from a structure that is equilibrated with the force field in the ensemble of interest, and that is representative of the reactive state you want to study. In this section we run a standard Amber protocol, check the trajectory and extract a snapshot with the solute imaged in the center of the box, in the format and with the periodic cell that CP2K expects.</p>

<br/>

<h3> <p align="center"><b> Equilibration protocol </b></p> </h3>

<br/>

| Stage | Input | Ensemble | Length | Restraints on solute heavy atoms |
| :--- | :--- | :--- | :--- | :--- |
| Minimization 1 | min1.in | - | 5000 steps | 10 kcal mol<sup>-1</sup> &Aring;<sup>-2</sup> |
| Minimization 2 | min2.in | - | 10000 steps | none |
| Heating 0 -> 300 K | heat.in | NVT, Langevin | 100 ps | 5 kcal mol<sup>-1</sup> &Aring;<sup>-2</sup> |
| Density equilibration | npt.in | NPT, Langevin + Monte Carlo barostat | 1 ns | 1 kcal mol<sup>-1</sup> &Aring;<sup>-2</sup> |
| Production | prod.in | NPT | 100 ns | none |

<br/>

- <p align="justify">The solute is defined with the mask !(:WAT,Na+,Cl-) &amp; !@H=, i.e. all heavy atoms of the protein, ligand and anything else that is not water or bulk ions. The MD stages use a 2 fs time step with SHAKE on the bonds to hydrogen (ntc=2, ntf=2), a 10 &Aring; cutoff and PME.</p>

<pre style="color: white; background-color: black;">
Heating 0 -> 300 K in NVT, solute heavy atoms restrained
&cntrl
  imin=0, irest=0, ntx=1,
  nstlim=50000, dt=0.002,
  ntc=2, ntf=2, cut=10.0, ntb=1,
  ntt=3, gamma_ln=2.0, ig=-1,
  tempi=0.0, temp0=300.0,
  nmropt=1,
  ntr=1, restraint_wt=5.0,
  restraintmask='!(:WAT,Na+,Cl-) & !@H=',
  ntpr=500, ntwx=5000, ntwr=5000,
  iwrap=1,
/
&wt type='TEMP0', istep1=0, istep2=40000, value1=0.0, value2=300.0, /
&wt type='TEMP0', istep1=40001, istep2=50000, value1=300.0, value2=300.0, /
&wt type='END', /
</pre>

<br/>

- <p align="justify">The run_md.sh script runs the stages in sequence, passing the restart file of each stage to the next one and using it as reference for the positional restraints. Use pmemd.cuda when a GPU is available (the default), or pmemd.MPI/sander with -n for the number of MPI tasks. The -l option overrides the number of steps of every stage, which is useful to test the whole chain in a few minutes before submitting the real job.</p>

<pre style="color: white; background-color: black;">
./run_md.sh -h
./run_md.sh -p complex.prmtop -c complex.rst7 -e pmemd.cuda
./run_md.sh -p complex.prmtop -c complex.rst7 -e sander -s "min1 min2 heat npt" -l 300   # quick test
</pre>

<pre style="color: white; background-color: black;">
[2026-10-25 10:55:06] running min1 with sander (input: complex.rst7)
[2026-10-25 10:56:32] running min2 with sander (input: min1.rst7)
[2026-10-25 10:58:04] running heat with sander (input: min2.rst7)
...
</pre>

<br/>

- <p align="justify">Production length depends on the question. For a QM/MM study of a reaction, you want enough sampling to know which conformations of the active site are populated and whether the reactive arrangement (the "near-attack conformation") is one of them. Running replicas with different initial velocities (ig=-1 gives a random seed) is usually more informative than a single long trajectory.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Checking the trajectory </b></p> </h3>

<br/>

- <p align="justify">Before choosing structures, check the energy, temperature, pressure and density of the equilibration stages (process_mdout.perl, or grep the mdout files), and then the stability of the protein and of the active site. analysis.cpptraj computes the backbone RMSD, the RMSD of the ligand without refitting (it measures how much it moves inside the protein) and the distances that define the catalytic triad and the position of the ligand relative to Ser195:</p>

<pre style="color: white; background-color: black;">
parm complex.prmtop
trajin prod.nc
autoimage anchor :1-291
rms bb_rmsd :1-290@CA,C,N first out rmsd.dat
rms lig_rmsd :291&amp;!@H= first nofit out rmsd.dat
distance d_ser_his  :206@OG  :44@NE2  out distances.dat
distance d_his_asp  :44@ND1  :100@OD1,OD2 out distances.dat
distance d_ser_lig  :206@OG  :291@O2  out distances.dat
run
</pre>

<pre style="color: white; background-color: black;">
cpptraj -i analysis.cpptraj
</pre>

<br/>

- <p align="justify">Choose the distances that matter for your mechanism: the nucleophile-electrophile distance and the angle of attack, the proton transfer distances, the coordination of metal ions, the hydrogen bonds of the oxyanion hole, etc. The histograms of these distances tell you how often the active site is ready for the reaction.</p>

<br/>

---

<br/>

<h3> <p align="center"><b> Selecting and imaging a snapshot </b></p> </h3>

<br/>

- <p align="justify">A common choice is the representative structure of the most populated cluster of the active site conformations (cluster.cpptraj), or several structures from different clusters if you want to estimate the conformational dependence of the barrier, or a frame that satisfies geometric criteria for the reaction.</p>

<pre style="color: white; background-color: black;">
parm complex.prmtop
trajin prod.nc 1 last 10
autoimage anchor :1-291
cluster c1 kmeans clusters 5 randompoint maxit 500 \
    rms :44,100,206,291&amp;!@H= \
    out cluster_vs_time.dat summary cluster_summary.dat info cluster_info.dat \
    repout cluster_rep repfmt restart
run
</pre>

<pre style="color: white; background-color: black;">
cat cluster_summary.dat
#Cluster   Frames     Frac  AvgDist    Stdev Centroid AvgCDist
       0        4    0.400    0.118    0.026        6    0.173
       ...
</pre>

<br/>

- <p align="justify">Two details are important for CP2K. First, the solute has to be whole and in the center of the box: autoimage puts the anchor molecule (the protein and the ligand, :1-291) in the center and wraps the solvent around it, so that the QM region is not split by the periodic boundary. Second, do not superimpose the frames on a reference before writing them (a "rms ... first" action without nofit): the rotation would misalign the system with the periodic box. The rms metric of the cluster command superimposes the frames internally and leaves the coordinates untouched.</p>

- <p align="justify">When you already know which frame you want, extract_snapshot.cpptraj writes it directly. The rst7 format keeps the box dimensions in the last line, and they are the ones to use in the &CELL section of CP2K (in NPT, the box of each frame is different from the initial one):</p>

<pre style="color: white; background-color: black;">
parm complex.prmtop
trajin prod.nc 2500 2500
autoimage anchor :1-291
trajout snapshot.rst7 restart
trajout snapshot.pdb pdb
run
</pre>

<pre style="color: white; background-color: black;">
tail -n 1 snapshot.rst7
  77.6793442  82.5830154  81.2552414  90.0000000  90.0000000  90.0000000
</pre>

<br/>

- <p align="justify">From now on we call the selected snapshot AMBER.rst7, the name used in the CP2K inputs of this tutorial and of the TS tutorial. The matching topology, AMBER.prmtop, is written in Section IV from complex.prmtop, after renaming the atom types that CP2K cannot handle. Open the snapshot with the topology in VMD and check that the protein and the ligand are whole and centered.</p>

<pre style="color: white; background-color: black;">
cp cluster_rep.c0.rst7 AMBER.rst7
</pre>

<br/>

<p align="justify">The input files of this section are in <a href="AMBER_INPUTS.zip">AMBER_INPUTS.zip</a>.</p>

<br/>
