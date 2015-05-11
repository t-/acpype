# Testing ACPYPE _amb2gmx_ function #

## Summarising ##

  * Using ACPYPE to convert AMBER _prmtop_ and _inpcrd_ files to GROMACS and the relative error for the total potential energy for a single point calculation is inferior to 0.009%, at the same range of error of [ffAMBER](http://ffamber.cnsm.csulb.edu/#validation).
  * Using **amb2gmx.pl** or even the modified **amb2gmx\_dihe.pl** (kindly provided by Guillem Portella), and the relative error is at least 0.017%, roughly twofold above ACPYPE's error.
  * If converting AMBER improper and proper dihedrals to GROMACS, the latter by using the multiple terms instead of RBs, and using `funct = 1` for all improper and proper dihedrals (instead of `9` and `4`), the results are simply the same.
  * ACPYPE is an alternative way to use PARMBSC0 in GROMACS.


**Table Sander ff99SB x ff99BSC0**
| **kCal/mol** | **Sander 99SB (_R1_)** | **Sander 99BSC0 (_R2_)** | **|ERROR|%** |
|:-------------|:-----------------------|:-------------------------|:-------------|
| **BOND**     |             94.3524|               94.3524|           0|
| **ANGLE**    |            829.1171|              829.1171|           0|
| **DIHEDRAL** |             462.101|              468.5221| 1.370500986|
| **VDW14**    |            263.4678|             2 63.4678|           0|
| **QQ14**     |          -2981.0039|            -2981.0039|           0|
| **VDW**      |           -366.7021|             -366.7021|           0|
| **QQ**       |           3650.5729|             3650.5729|           0|
| **TOT POT**  |           1951.9052|             1958.3262| 0.327882045|


**Table ff99SB**
| **kJ/mol**   | **Sander (_R1_)** | **ACPYPE (_R3_)** | **|ERROR|%**   | **AMB2GMX (_R5_)** | **|ERROR|%**    | **GMX45 (_R7_)** | **|ERROR|%**   |
|:-------------|:------------------|:------------------|:---------------|:-------------------|:----------------|:-----------------|:---------------|
| **BOND**     | 394.7704416   | 394.773       | 0.000648069  | 394.773        | 0.000648069   | 394.772      | 0.00039476   |
| **ANGLE**    | 3469.025946   | 3469.02       | 0.000171414  | 3469.02        | 0.000171414   | 3469.02      | 0.000171414  |
| **DIHEDRAL** | 1933.430584   | 1933.42546    | 0.000265021  | 1934.24        | **0.04184672**  | 1933.4321    | 7.84098E-05  |
| **VDW14**    | 1102.349275   | 1102.35       | 6.57504E-05  | 1102.35        | 6.57504E-05   | 1102.35      | 6.57504E-05  |
| **QQ14**     | -12472.52032  | -12472.5      | -0.000162899 | -12472.4       | -0.000964671  | -12472.5     | -0.000162899 |
| **VDW**      | -1534.281586  | -1534.28      | -0.000103397 | -1534.28       | -0.000103397  | -1534.28     | -0.000103397 |
| **QQ**       | 15273.99701   | 15274.7       | 0.004602293  | 15274.5        | 0.003292981   | 15274.5      | 0.003292981  |
| **TOT POT**  | 8166.771357   | 8167.49       | 0.008798826  | 8168.17        | **0.017123091** | 8167.36      | 0.007207264  |


**Table ff99BSC0**
| **kJ/mol**   | **Sander (_R2_)** | **ACPYPE (_R4_)** | **|ERROR|%**   | **AMB2GMX (_R6_)** | **|ERROR|%**    | **AMB2GMX\_DIHE (_R8_)** | **|ERROR|%**    | **ACPYPE GMX45 (_R9_)** | **|ERROR|%**   |
|:-------------|:------------------|:------------------|:---------------|:-------------------|:----------------|:-------------------------|:----------------|:------------------------|:---------------|
| **BOND**     | 394.7704416     | 394.773         | 0.000648069  | 394.773          | 0.000648069   | 394.773               | 0.000648069   | 394.773               | 0.000648069  |
| **ANGLE**    | 3469.025946     | 3469.02         | 0.000171414  | 3469.02          | 0.000171414   | 3469.02               | 0.000171414   | 3469.02               | 0.000171414  |
| **DIHEDRAL** | 1960.296466     | 1960.29         | 0.000329868  | 2198.1           | **10.81859486** | 1961.773              | **0.075265263** | 1960.29               | 0.000329868  |
| **VDW14**    | 1102.349275     | 1102.35         | 6.57504E-05  | 1102.35          | 6.57504E-05   | 1102.35               | 6.57504E-05   | 1102.35               | 6.57504E-05  |
| **QQ14**     | -12472.52032    | -12472.5        | -0.000162899 | -12472.4         | -0.000964671  | -12472.4              | -0.000964671  | -12472.5              | -0.000162899 |
| **VDW**      | -1534.281586    | -1534.28        | -0.000103397 | -1534.28         | -0.000103397  | -1534.28              | -0.000103397  | -1534.28              | -0.000103397 |
| **QQ**       | 15273.99701     | 15274.7         | 0.004602293  | 15274.5          | 0.003292981   | 15274.5               | 0.003292981   | 15274.7               | 0.004602293  |
| **TOT POT**  | 8193.636821     | 8194.35         | 0.008703304  | 8432.04          | **2.82734877**  | 8195.71               | **0.025295907** | 8194.35               | 0.008703304  |


## Detailing ##

Please see the script and its comments below. You may be able to reproduce the results below as long as you have the programmes needed.

```
# $AMBERHOME/dat/leap/parm/frcmod.parmbsc0: fixed CI mass to 12.01
# def error(a,b): return abs(a - b)/max(a,b) *100
# using gromacs 4.5; amb2gmx_dihe.pl provided by Guillem Portella

wget -c "http://www.pdbe.org/download/1BNA" -O 1BNA.pdb

# Assuming DNA.pdb (= 1BNQ.pdb)
grep 'ATOM  ' 1BNA.pdb >| DNA.pdb

# Create a input file with commands for tleap and run it
cat << EOF >| leap.in
source leaprc.ff99SB
dna = loadpdb DNA.pdb
saveamberparm dna DnaAmberSB.prmtop DnaAmber.inpcrd
source leaprc.ff99bsc0
dna = loadpdb DNA.pdb
saveamberparm dna DnaAmberBSC0.prmtop DnaAmber.inpcrd
savepdb dna DnaAmber.pdb
quit
EOF
tleap -f leap.in >| leap.out

cat << EOF >| mdin
Single point
&cntrl
imin=0, maxcyc=0,
ntmin=2,
ntb=0,
igb=0,
cut=999,
//
EOF

# Round off DnaAmber.inpcrd to 3 decimals to equals DnaAmber.pdb
python
inp = open('DnaAmber.inpcrd').readlines()
new = inp[:2]
for line in inp[2:]:
    d =["%5.7f" % x for x in [round(x,3) for x in list(map(float, [line[i:i + 12] for i in range(0, 72, 12)]))]]
    new.append("%12s" * 6 % tuple(d) + '\n')

open('DnaAmber.inpcrd','w').writelines(new)
exit()

###################
### SANDER 99SB ###
###################
sander -O -i mdin -o mdout -p DnaAmberSB.prmtop -c DnaAmber.inpcrd
\mv mdinfo DnaAmberSB.mdinfo
cat DnaAmberSB.mdinfo
# (R1)
 NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
 Etot   =      1951.9052  EKtot   =         0.0000  EPtot      =      1951.9052
 BOND   =        94.3524  ANGLE   =       829.1171  DIHED      =       462.1010
 1-4 NB =       263.4678  1-4 EEL =     -2981.0039  VDWAALS    =      -366.7021
 EELEC  =      3650.5729  EHBOND  =         0.0000  RESTRAINT  =         0.0000

#####################
### SANDER 99BSC0 ###
#####################
sander -O -i mdin -o mdout -p DnaAmberBSC0.prmtop -c DnaAmber.inpcrd
\mv mdinfo DnaAmberBSC0.mdinfo
cat DnaAmberBSC0.mdinfo
# (R2)
 NSTEP =        0   TIME(PS) =       0.000  TEMP(K) =     0.00  PRESS =     0.0
 Etot   =      1958.3262  EKtot   =         0.0000  EPtot      =      1958.3262
 BOND   =        94.3524  ANGLE   =       829.1171  DIHED      =       468.5221
 1-4 NB =       263.4678  1-4 EEL =     -2981.0039  VDWAALS    =      -366.7021
 EELEC  =      3650.5729  EHBOND  =         0.0000  RESTRAINT  =         0.0000

### NOTE: R1 x R2: only DIHED diffs (ERRdih = |1.371| %, ERRpot = |0.328| % )

# Create SPE.mdp file #single point energy
cat << EOF >| SPE.mdp
define = -DFLEXIBLE
integrator               = md
nsteps                   = 0
dt                       = 0.001
constraints              = none
emtol                    = 10.0
emstep                   = 0.01
nstcomm                  = 1
ns_type                  = simple
nstlist                  = 0
rlist                    = 0
rcoulomb                 = 0
rvdw                     = 0
Tcoupl                   = no
Pcoupl                   = no
gen_vel                  = no
nstxout                  = 1
pbc                      = no
nstlog = 1
nstenergy = 1
nstvout = 1
nstfout = 1
nstxtcout = 1
comm_mode = ANGULAR
continuation = yes
EOF

###################
### ACPYPE 99SB ### DIHE FUNCT 3 (amber prop. dihe -> RB), 1 (amber impr. dihe.)
###################
acpype -p DnaAmberSB.prmtop -x DnaAmber.inpcrd
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmberSB_GMX.top -o speSB.tpr
mdrun -v -deffnm speSB
echo 1 2 3 4 5 6 7 8 9 | g_energy -f speSB.edr
# (R3)
Bond                        394.773         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Proper Dih.                 2.64546         --          0          0  (kJ/mol)
Ryckaert-Bell.              1930.78         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.5         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.7         --          0          0  (kJ/mol)
Potential                   8167.49         --          0          0  (kJ/mol)

### NOTE: R3 x R1: DnaAmberSB_GMX x DnaAmberSB.mdinfo |errors|: bond=0.001, ang=0.000, dih=0.000, vdW14=0.000, QQ14=0.000, vdW=0.000, QQ=0.005, pot=0.009

#####################
### ACPYPE 99BSC0 ### DIHE DIHE FUNCT 3 (amber prop. dihe. -> RB), 1 (amber impr. dihe. and bsc0 dihe.)
#####################
acpype -p DnaAmberBSC0.prmtop -x DnaAmber.inpcrd
\mv DnaAmberBSC0_GMX.top DnaAmberBSC0_GMX_31.top
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmberBSC0_GMX_31.top -o speBSC0_31.tpr
mdrun -v -deffnm speBSC0_31
echo 1 2 3 4 5 6 7 8 9 | g_energy -f speBSC0_31.edr
# (R4)
Bond                        394.773         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Proper Dih.                  135.56         --          0          0  (kJ/mol)
Ryckaert-Bell.              1824.73         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.5         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.7         --          0          0  (kJ/mol)
Potential                   8194.35         --          0          0  (kJ/mol)

### NOTE: R4 x R2: DnaAmberBSC0_GMX x DnaAmberBSC0.mdinfo |errors|: bond=0.001, ang=0.000, dih=0.000, vdW14=0.000, QQ14=0.000, vdW=0.000, QQ=0.005, pot=0.009 (same values as NOTE: R3 x R1)

####################
### AMB2GMX 99SB ###  DIHE DIHE FUNCT 3, 3 (RB for all dihe.)
####################
~/Programmes/amb2gmx.pl --prmtop DnaAmberSB.prmtop --crd DnaAmber.inpcrd --outname DnaAmb2gmxSB --debug
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmb2gmxSB.top -o speAmb2gmxSB.tpr
mdrun -v -deffnm speAmb2gmxSB
echo 1 2 3 4 5 6 7 8 | g_energy -f speAmb2gmxSB.edr
# (R5)
Bond                        394.773         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Ryckaert-Bell.              1934.24         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.4         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.5         --          0          0  (kJ/mol)
Potential                   8168.17         --          0          0  (kJ/mol)

### NOTE: R5 x R1: DnaAmb2gmxSB x DnaAmberSB.mdinfo |errors| : bond=0.001, ang=0.000, dih=0.042, vdW14=0.000, QQ14=0.001, vdW=0.000, QQ=0.003, pot=0.017 (only dih, QQ14, QQ diffs from NOTE: R3 x R1)

######################
### AMB2GMX 99BSC0 ### DIHE FUNCT 3, 3 (RB for all dihe.)
######################
~/Programmes/amb2gmx.pl --prmtop DnaAmberBSC0.prmtop --crd DnaAmber.inpcrd --outname DnaAmb2gmxBSC0 --debug
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmb2gmxBSC0.top -o speAmb2gmxBSC0.tpr
mdrun -v -deffnm speAmb2gmxBSC0
echo 1 2 3 4 5 6 7 8 | g_energy -f speAmb2gmxBSC0.edr
# (R6)
Bond                        394.773         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Ryckaert-Bell.               2198.1         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.4         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.5         --          0          0  (kJ/mol)
Potential                   8432.04         --          0          0  (kJ/mol)

### NOTE: R6 x R2: DnaAmb2gmxBSC0 x DnaAmberBSC0.mdinfo |errors|: bond=0.001, ang=0.000, dih=10.819, vdW14=0.000, QQ14=0.001, vdW=0.000, QQ=0.003, pot=2.827 (only dih, QQ14, QQ diffs from NOTE: R3 x R1)

##################
### GMX45 99SB ### DIHE FUNCT 9 (amber prop. dihe.), 4 (amber impr. dihe.)
##################
pdb2gmx -f DnaAmber.pdb -o DnaAmberSBGMX45.pdb -ff amber99sb -water none -p DnaAmberSBGMX45
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmberSBGMX45.top -o speSBGMX45.tpr
mdrun -v -deffnm speSBGMX45
echo 1 2 3 4 5 6 7 8 9 | g_energy -f speSBGMX45.edr
# (R7)
Bond                        394.772         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Proper Dih.                 1930.78         --          0          0  (kJ/mol)
Improper Dih.                2.6521         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.5         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.5         --          0          0  (kJ/mol)
Potential                   8167.36         --          0          0  (kJ/mol)

### NOTE: R7 x R1: DnaAmberSBGMX45 x DnaAmberSB.mdinfo |errors|: bond=0.000, ang=0.000, dih=0.000, vdW14=0.000, QQ14=0.000, vdW=0.000, QQ=0.003, pot=0.007 (only bond and QQ diffs from NOTE: R3 x R1)

###########################
### AMB2GMX_DIHE 99BSC0 ### DIHE FUNCT 3 (even for amber impr. dihe.), 1 (for bsc0 dihe.)
###########################
~/Programmes/amb2gmx_dihe.pl --prmtop DnaAmberBSC0.prmtop --crd DnaAmber.inpcrd --outname DnaAmb2gmx_diheBSC0 --debug
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmb2gmx_diheBSC0.top -o DnaAmb2gmx_diheBSC0.tpr
mdrun -v -deffnm DnaAmb2gmx_diheBSC0
echo 1 2 3 4 5 6 7 8 9 | g_energy -f DnaAmb2gmx_diheBSC0.edr
# (R8)
Bond                        394.773         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Proper Dih.                 133.593         --          0          0  (kJ/mol)
Ryckaert-Bell.              1828.18         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.4         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.5         --          0          0  (kJ/mol)
Potential                   8195.71         --          0          0  (kJ/mol)

### NOTE: R8 x R2: DnaAmb2gmx_diheBSC0 x DnaAmberBSC0.mdinfo |errors|: bond=0.001, ang=0.000, dih=0.075, vdW14=0.000, QQ14=0.001, vdW=0.000, QQ=0.003, pot=0.025 (only dihe QQ14 and QQ diffs from NOTE: R3 x R1)

###########################
### ACPYPE GMX45 99BSC0 ### DIHE FUNCT 9 (amber prop. dihe.), 4 (amber impr. dihe. and bsc0 dihe.)
###########################
acpype -p DnaAmberBSC0.prmtop -x DnaAmber.inpcrd -r
\mv DnaAmberBSC0_GMX.top DnaAmberBSC0_GMX_94.top
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmberBSC0_GMX_94.top -o speBSC0_94.tpr
mdrun -v -deffnm speBSC0_94
echo 1 2 3 4 5 6 7 8 9 | g_energy -f speBSC0_94.edr
# (R9)
Bond                        394.773         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Proper Dih.                 1824.73         --          0          0  (kJ/mol)
Improper Dih.                135.56         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.5         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.7         --          0          0  (kJ/mol)
Potential                   8194.35         --          0          0  (kJ/mol)

# Replacing in DnaAmberBSC0_GMX.top dihedrals funct 9 and 4 by 1 and 1 (DIHE FUNCT 1, 1)
grompp -f SPE.mdp -c DnaAmber.pdb -p DnaAmberBSC0_GMX_11.top -o speBSC0_11.tpr
mdrun -v -deffnm speBSC0_11
echo 1 2 3 4 5 6 7 8 | g_energy -f speBSC0_11.edr
# (R10)
Bond                        394.773         --          0          0  (kJ/mol)
Angle                       3469.02         --          0          0  (kJ/mol)
Proper Dih.                 1960.29         --          0          0  (kJ/mol)
LJ-14                       1102.35         --          0          0  (kJ/mol)
Coulomb-14                 -12472.5         --          0          0  (kJ/mol)
LJ (SR)                    -1534.28         --          0          0  (kJ/mol)
Coulomb (SR)                15274.7         --          0          0  (kJ/mol)
Potential                   8194.35         --          0          0  (kJ/mol)

### NOTE: R4, R9 and R10 give absolutely the same results!!!
### R10 also tested with Gromacs 4.0.5 and gave the same result
```