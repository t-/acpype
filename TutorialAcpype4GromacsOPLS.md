# Tutorial Using ACPYPE for GROMACS with OPLS/AA forcefield #

## Introduction ##

This tutorial is to show how to prepare a system to run on GROMACS with OPLS/AA
forcefield, starting with a PDB file for a complex protein/ligand.

It is a mere proof of concept and OPLS/AA is very experimental in ACPYPE.
One should see how ACPYPE is trying to do it in HowAcpypeWorks.
If you have suggestions about how to improve this tutorial, please send a
comment (at the bottom of the page).

**NB:** Besides **acpype**, **antechamber** and **babel**, you will need GROMACS.

## Getting GROMACS ##

Install [GROMACS](http://www.gromacs.org/). The current version is 4.0.7.
Something like:

  * `sudo apt-get install gromacs` # if you use Ubuntu Linux, or

  * `fink install gromacs` # if you use Mac

should do the trick.

## Running an Example ##

This is for protein 1BVG.pdb (get it at [PDB](http://www.pdb.org)), a homodimer
(HIV protease) with a ligand called DMP. We will use force field OPLS/AA.

Luckily, this pdb file has all hydrogens for the ligand, which is necessary for
**antechamber**. One can use either, e.g., `babel -h _mol_w/o_H_.pdb _mol_with_H.pdb`
or [YASARA View](http://www.yasara.org) to automatically add missing hydrogens to
your compound. The former just puts 'H' for atom names while the latter puts
more meaningful atom name, e.g., 'HCA' for a H bonded to a CA and not a simply
'H' as **babel** does.

In a script-like way:
```
# Assuming Complex.pdb (= 1BVG.pdb), split it in Protein.pdb and Ligand.pdb
wget -c "http://www.pdb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=1BVG" -O 1BVG.pdb
grep 'ATOM  ' 1BVG.pdb>| Protein.pdb
grep 'HETATM' 1BVG.pdb>| Ligand.pdb

# Process with pdb2gmx and define water
pdb2gmx -ff oplsaa -f Protein.pdb -o Protein2.pdb -p Protein.top -water spce -ignh

# Generate Ligand topology file with acpype (GAFF)
acpype -i Ligand.pdb

# Merge Protein2.pdb + updated Ligand_NEW.pdb -> Complex.pdb
grep -h ATOM Protein2.pdb Ligand.acpype/Ligand_NEW.pdb >| Complex.pdb

# Edit Protein.top -> Complex.top
\cp Ligand.acpype/Ligand_GMX_OPLS.itp Ligand.itp
\cp Protein.top Complex.top
# See NB(1) below
cat Complex.top | sed '/\;\ Include\ chain\ topologies/i\
#include "Ligand.itp"
' >| Complex2.top
echo "Ligand   1" >> Complex2.top
\mv Complex2.top Complex.top

# Setup the box and add water
editconf -bt triclinic -f Complex.pdb -o Complex.pdb -d 1.0
genbox -cp Complex.pdb -cs -o Complex_b4ion.pdb -p Complex.top

# Create em.mdp file
cat << EOF >| em.mdp
define                   = -DFLEXIBLE
integrator               = cg ; steep
nsteps                   = 200
constraints              = none
emtol                    = 1000.0
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

# Need to comment out some dih parameters in Ligand.itp that OPLS/AA does not recognise.
# Note that the parameters values here are derived for AMBER99SB/GAFF and NOT for OPLS/AA.
# A proper solution would involve finding the right OPLS atom types and parameters
# calculation according to [http://dx.doi.org/10.1021%2Fja9621760 Jorgensen et al. (1996)].
# Be aware of it!

# To illustrate this problem, try the 1st command in 'Setup ions' session without
# the fix right below.

sed s/15\ \ 14\ \ \ 2\ \ 36\ \ \ 3\ \;/15\ \ 14\ \ \ 2\ \ 36\ \ \ 3\ \ / Ligand.itp |\
sed s/22\ \ 21\ \ \ 5\ \ 37\ \ \ 3\ \;/22\ \ 21\ \ \ 5\ \ 37\ \ \ 3\ \ / | \
sed s/36\ \ \ 1\ \ 37\ \ \ 5\ \ \ 3\ \;/36\ \ \ 1\ \ 37\ \ \ 5\ \ \ 3\ \ / | \
sed s/36\ \ \ 1\ \ 37\ \ 28\ \ \ 3\ \;/36\ \ \ 1\ \ 37\ \ 28\ \ \ 3\ \ / | \
sed s/36\ \ \ 2\ \ \ 3\ \ 39\ \ \ 3\ \;/36\ \ \ 2\ \ \ 3\ \ 39\ \ \ 3\ \ / | \
sed s/37\ \ \ 1\ \ 36\ \ \ 2\ \ \ 3\ \;/37\ \ \ 1\ \ 36\ \ \ 2\ \ \ 3\ \ / | \
sed s/37\ \ \ 1\ \ 36\ \ \ 6\ \ \ 3\ \;/37\ \ \ 1\ \ 36\ \ \ 6\ \ \ 3\ \ / | \
sed s/37\ \ \ 5\ \ \ 4\ \ 40\ \ \ 3\ \;/37\ \ \ 5\ \ \ 4\ \ 40\ \ \ 3\ \ / | \
sed s/\;\ \ \ 180\.00/\ \ \ \ 180\.00/g > TemPfiLe; \mv -b TemPfiLe Ligand.itp


# Setup ions
grompp -f em.mdp -c Complex_b4ion.pdb -p Complex.top -o Complex_b4ion.tpr
\cp Complex.top Complex_ion.top
echo 13| genion -s Complex_b4ion.tpr -o Complex_b4em.pdb -neutral -conc 0.15 -p Complex_ion.top -norandom -nname CL- -pname NA+
\mv Complex_ion.top Complex.top

# Run minimisaton
grompp -f em.mdp -c Complex_b4em.pdb -p Complex.top -o em.tpr
mdrun -v -deffnm em

# Run a short simulation
grompp -f md.mdp -c em.gro -p Complex.top -o md.tpr
mdrun -v -deffnm md

# or with openmpi, for a dual core
grompp -f em.mdp -c Complex_b4em.pdb -p Complex.top -o em.tpr
om-mpirun -n 2 mdrun_mpi -v -deffnm em
grompp -f md.mdp -c em.gro -p Complex.top -o md.tpr
om-mpirun -n 2 mdrun_mpi -v -deffnm md

# Visualise with VMD
vmd md.gro md.trr
```

**NB(1):** `#include "Ligand.itp"` has to be inserted right after `ffamber**.itp`
line and before `Protein_*.itp` line in _Complex.top_.

Voila!