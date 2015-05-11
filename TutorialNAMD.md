# How To Use NAMD with AMBER Output Topologies From ACPYPE #

## Introduction ##

This tutorial is to show how to prepare a system to run on NAMD, starting
with a PDB file for a complex protein/ligand.

It is a mere proof of concept.
If you have suggestions about how to improve this tutorial, please send a
comment (at the bottom of the page).

For this tutorial I use [NAMD](http://www.ks.uiuc.edu/Research/namd). It is free
and works fine with AMBER input files, however with some caveats:

  * Polarizable parameters in AMBER are not supported.
  * NAMD does not support the 10-12 potential terms in some old AMBER versions. When non-zero 10-12 parameter is encountered in PARM file, NAMD will terminate.
  * NAMD has several exclusion policy options, defined by `exclude`. The way AMBER dealing with exclusions corresponds to the "scaled1-4" in NAMD. So for simulations using AMBER force field, one would specify "exclude scaled1-4" in the configuration file, and set `1-4scaling` to the inverse value of SCEE as would be used in AMBER.
  * NAMD does not read periodic box lengths in PARM or coordinate file. They must be explicitly specified in NAMD configuration file.
  * By default, NAMD applies switching functions to the non-bond interactions within the cutoff distance, which helps to improve energy conservation, while AMBER does not use switching functions so it simply truncates the interactions at cutoff. However, if "authentic" AMBER cutoff simulations are desired, the switching functions could be turned off by specifying "switching off" in NAMD configuration file.
  * NAMD and AMBER may have different default values for some parameters (e.g., the tolerance of SHAKE). One should check [NAMD manual](http://www.ks.uiuc.edu/Research/namd/current/ug/node14.html) for accurate descriptions of the NAMD options.

**NB:** One doesn't really need **acpype**, only **antechamber** to do it, but I
believe **acpype** makes this chore easier.

## Running an Example ##

This is for protein 1BVG.pdb (get it at [PDB](http://www.pdb.org)), a homodimer
(HIV protease) with a ligand called DMP. We will use force field AMBER99SB.

Luckily, this pdb file has all hydrogens for the ligand, which is necessary for
**antechamber**. One can use either, e.g., `babel -h _mol_w/o_H_.pdb _mol_with_H.pdb`
or [YASARA View](http://www.yasara.org) to automatically add missing hydrogens to
your compound. The former just puts 'H' for atom names while the latter puts
more meaningful atom name, e.g., 'HCA' for a H bonded to a CA and not a simply
'H' as **babel** does.

In a script-like way:
```
# Assuming 1BVG.pdb, create Ligand.pdb
wget -c "http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=1BVG" -O 1BVG.pdb
grep 'HETATM' 1BVG.pdb>| Ligand.pdb

# Edit 1BVG.pdb and create ComplexAmber.pdb
sed s/HIS\ /HID\ /g 1BVG.pdb | sed s/H1\ \ PRO/H3\ \ PRO/g >| ComplexAmber.pdb

# Generate Ligand topology file with acpype (GAFF)
acpype -i Ligand.pdb

# Create a input file with commands for tleap and run it
cat << EOF >| leap.in
verbosity 1
source leaprc.ff99SB
source leaprc.gaff
loadoff Ligand.acpype/Ligand_AC.lib
loadamberparams Ligand.acpype/Ligand_AC.frcmod
complex = loadpdb ComplexAmber.pdb
solvatebox complex TIP3PBOX 10.0
addions complex Na+ 23
addions complex Cl- 27
saveamberparm complex ComplexAmber.prmtop ComplexAmber.inpcrd
savepdb complex ComplexNAMD.pdb
quit
EOF
tleap -f leap.in >| leap.out

# Create NAMD run file
cat << EOF >| run_namd.conf
rigidBonds     all
rigidTolerance 0.0005  # Default is  0.00001
outputEnergies 50  # Energy output frequency
DCDfreq        2  # Trajectory file frequency
timestep       2  # in unit of fs
temperature    300  # Initial temp for velocity assignment
cutoff         10
switching      off  # Turn off the switching functions
PME            on  # Use PME for electrostatic calculation
# Orthogonal periodic box size
cellBasisVector1   81.8922960  0  0
cellBasisVector2   0  63.8569220  0
cellBasisVector3   0  0  60.8145690
PMEGridSizeX   82
PMEGridSizeY   64
PMEGridSizeZ   61
# NAMD doesn't force neutralisation of charge
amber          on  # Specify this is AMBER force field
parmfile       ComplexAmber.prmtop  # Input PARM file
ambercoor      ComplexAmber.inpcrd  # Input coordinate file
outputname     run_namd  # Prefix of output files
exclude        scaled1-4
1-4scaling     0.833333  # =1/1.2, default is 1.0
minimize       200
reinitvels     300
run            1000 ;# 2ps
EOF

# Run NAMD
namd2 run_namd.conf >| run_namd.log

# if dual core is available:
namd2 +p2 run_namd.conf >| run_namd2.log

# Visualise with VMD
vmd -parm7 ComplexAmber.prmtop -rst7 ComplexAmber.inpcrd run_namd.dcd
```

## Comparing NAMD with GROMACS ##

Now, let's take AMBER input files (`ComplexAmber.prmtop` and `ComplexAmber.inpcrd`)
and convert them to GROMACS input files right for simulations for the very same
molecular system.

In a script-like way:
```
acpype -p ComplexAmber.prmtop -x ComplexAmber.inpcrd

# Create em.mdp file
cat << EOF >| em.mdp
define                   = -DFLEXIBLE
integrator               = cg ; steep
nsteps                   = 200
constraints              = none
emtol                    = 10.0
nstcgsteep               = 10 ; do a steep every 10 steps of cg
emstep                   = 0.01 ; used with steep
nstcomm                  = 1
coulombtype              = PME
ns_type                  = grid
rlist                    = 1.0
rcoulomb                 = 1.0
rvdw                     = 1.4
Tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
nstxout                  = 0 ; write coords every # step
optimize_fft             = yes
EOF

# Create md.mdp file
cat << EOF >| md.mdp
integrator               = md
nsteps                   = 1000
dt                       = 0.002
constraints              = all-bonds
nstcomm                  = 1
ns_type                  = grid
rlist                    = 1.2
rcoulomb                 = 1.1
rvdw                     = 1.0
vdwtype                  = shift
rvdw-switch              = 0.9
coulombtype              = PME-Switch
Tcoupl                   = v-rescale
tau_t                    = 0.1 0.1
tc-grps                  = protein non-protein
ref_t                    = 300 300
Pcoupl                   = parrinello-rahman
Pcoupltype               = isotropic
tau_p                    = 0.5
compressibility          = 4.5e-5
ref_p                    = 1.0
gen_vel                  = yes
nstxout                  = 2 ; write coords every # step
lincs-iter               = 2
DispCorr                 = EnerPres
optimize_fft             = yes
EOF

# Run minimisaton
grompp -f em.mdp -c ComplexAmber_GMX.gro -p ComplexAmber_GMX.top -o em.tpr
mdrun -v -deffnm em

# Run a short simulation
grompp -f md.mdp -c em.gro -p ComplexAmber_GMX.top -o md.tpr
mdrun -v -deffnm md

# or with openmpi, for a dual core
grompp -f em.mdp -c ComplexAmber_GMX.gro -p ComplexAmber_GMX.top -o em.tpr
om-mpirun -n 2 mdrun_mpi -v -deffnm em
grompp -f md.mdp -c em.gro -p ComplexAmber_GMX.top -o md.tpr
om-mpirun -n 2 mdrun_mpi -v -deffnm md

# Visualise with VMD
vmd md.gro md.trr
```

Voila!

#### NOTE ####
  * I am trying to set up here the same simulation conditions for NAMD 2.7 and GROMACS 4.0.4, althought NAMD runs in double precision while GROMACS in single precision.
  * GROMACS is much faster than NAMD in my Mac Intel Duo Core 2Ghz:
    * 1 proc.: 694%
    * 2 proc.: 624%
  * But NAMD scales better when using 2 proc.:
    * gain of 80% for NAMD
    * gain of 64% for GROMACS