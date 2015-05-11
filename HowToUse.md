# How To Use ACPYPE #

## Introduction ##

To run **acpype** with its all functionalities, you need **ANTECHAMBER** from package
[AmberTools](http://ambermd.org/#AmberTools) and
[Open Babel](http://openbabel.org/wiki/Main_Page) if your input files are of PDB
format.

However, if one wants **acpype** just to emulate **amb2gmx.pl**, one needs nothing
at all but **[Python](http://www.python.org)**.

At the moment, **acpype** is only available for download via **svn**:

  * `svn checkout http://acpype.googlecode.com/svn/trunk/ acpype`

Yet, if some reason you cannot use **svn**, one still can get **acpype** with:

  * `wget http://acpype.googlecode.com/svn/trunk/acpype.py`

But be aware that one may run in extra troubles and I am not willing to support
this way.

## To Test ##

At folder **acpype/test**, type:

  * `../acpype.py -i FFF.pdb`

It'll create a folder called **FFF.acpype**, and inside it one may find topology
files for GROMACS and CNS/XPLOR.

To get help and more information, type:

  * `../acpype.py -h`

## To Install ##

At folder **acpype**, type:

  * `ln -s $PWD/acpype.py /usr/local/bin/acpype`

And re-login or start another shell session.

## To Verify with GMX ##
```
cd FFF.acpype/
grompp -c FFF_GMX.gro -p FFF_GMX.top -f em.mdp -o em.tpr
mdrun -v -deffnm em
# And if you have VMD
vmd em.gro em.trr

# For MD, do:
grompp -c em.gro -p FFF_GMX.top -f md.mdp -o md.tpr
mdrun -v -deffnm md
vmd md.gro md.trr

# with openmpi, for a dual core
grompp -c FFF_GMX.gro -p FFF_GMX.top -f em.mdp -o em.tpr
om-mpirun -n 2 mdrun_mpi -v -deffnm em
grompp -c em.gro -p FFF_GMX.top -f md.mdp -o md.tpr
om-mpirun -n 2 mdrun_mpi -v -deffnm md
vmd md.gro md.trr
```

## To Emulate _[amb2gmx.pl](http://ffamber.cnsm.csulb.edu/tools.html)_ ##

For any given **prmtop** and **inpcrd** files (outputs from AMBER LEaP), type:

  * ` acpype -p FFF_AC.prmtop -x FFF_AC.inpcrd`

The output files `FFF_GMX.gro` and `FFF_GMX.top` will be generated at the same
folder of the input files.

## To Verify with CNS/XPLOR ##
At folder **FFF.acpype**, type:

  * `cns < FFF_CNS.inp`

## To Verify with NAMD ##

  * see [TutorialNAMD](TutorialNAMD.md)