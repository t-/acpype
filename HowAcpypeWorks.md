# How ACPYPE Works #

It is totally based on **ANTECHAMBER** and although it may work with old versions,
it is recommended to use the one that comes in
**[AmberTools](http://ambermd.org/AmberTools-get.html)** (version 12 now).

It basically does what [amb2gmx](http://ffamber.cnsm.csulb.edu/tools.html)
does but with some remarkable differences:

  * as in **amb2gmx.pl**, torsionals (proper and improper) are treated as Ryckaert-Bellemans in GROMACS to use combine multiple AMBER torsions per quartet. However **acpype** follows the approach used in porting AMBER force fields to GROMACS, which  means that improper dihedrals are separated from proper dihedrals and treated as proper dihedrals in GROMACS to use correct AMBER analytical function;

  * **acpype** knows about PARMBSC0 dihedrals parameters and correctly convert them to GROMACS (see [TestingAcpypeAmb2gmx](TestingAcpypeAmb2gmx.md));

  * no dependence at all of **rdparm** or **ambpdb**. Although now **rdparm** is included within **AmberTools** (before it was only available within **AMBER** package), **ambpdb** is not, which means that if one wants to use **amb2gmx.pl**, one needs **AMBER#** package (where # can be 7 or bigger);

  * no fixed size limit in number of atoms as for **amb2gmx.pl**, since usual compilation of **rdparm** limits up to 150k atoms and 50k residues. **acpype** is only limited by the computer's memory;

  * it also reads and converts correctly box (or octahedron) parameters in INPCRD to GRO file, when available, otherwise new box parameters will be calculated. Besides it knows if TIP3P or SPCE water and applies the correct parameters;

  * it generates topologies for CNS/XPLOR and a tentative OPLS/AA for GROMACS (experimental);

  * **acpype** tries to smooth the integration with GROMACS for immediate use (see [TutorialAcpype4Gromacs](TutorialAcpype4Gromacs.md) and [TutorialAcpype4GromacsOPLS](TutorialAcpype4GromacsOPLS.md)).

**HINT:** One can use [YASARA AutoSMILES Server](http://www.yasara.org/autosmilesserver.htm) to obtain probably better charges and insert them in a _mol2_ file. If given this _mol2_ file as input to **acpype** and with option `-c user`, **acpype** will take these charges defined in the _mol2_ input file.

## IMPORTANT ##

Because ACPYPE is based on Antechamber, it will as well inherit some of its limitations like not being possible to work with organic molecule with open valences; containing others atoms than C, N, O, S, P, H, F, Cl, Br and I; or covalently bonded to another molecule. If one wants parameters for a modified amino acid residue, one way of getting it is by neutralising the N- and C- termini and then fit manually the additional parameters to the modified residue.

Ideally, one should know well about the molecule one is trying to parametrise
and should read the pertinent papers related to the force field intended to be
used (see references below). Therefore be **_warned_** if using
**acpype/antechamber** as a **_blackbox_**.

## A Comparative Test for ACPYPE ##

Running a testing routine, I used `pymol` to generate 22 pdb files.
They are tripeptides of every usual amino acid residue
including 2 variants for HIS. All peptides have N- and T- termini.

I used GROMACS 4.5 with all these tripeptides to generate topology files with AMBER99SB forcefield as reference.

OBS: peptides JJJ.pdb (Hip-Hip-Hip) and RRR.pdb (Arg-Arg-Arg): their net charge should be +3, but _gasteiger_
method (the way how ACPYPE try to guess the net charge of the molecule if not
given) failed to get the right value, instead it says it is 'Zero', for which SQM programme will fail to calculate the atomic partial charges. So be aware
that ACPYPE may guess a wrong net charge when running it on your own molecule.

Next, I used ACPYPE to generate the topologies and parameters for these 22 entries and compare theirs
results with the ones from GROMACS reference.

By using ACPYPE with option '-a amber' (which means parm99.dat + gaff.dat +
frcmod.ff99SB + frcmod.parmbsc0), it seems to give better results than using just GAFF (the default option in ACPYPE) when comparing with GROMACS outputs, i.e., the former basically got almost all atom types and parameters identical to the GROMACS with AMBER99SB reference, while the latter missed just a few more.

A more detailed comparison between ACPYPE with option '-a amber' against
GROMACS' AMBER99SB (including a Single Point Energy minimisation with GROMACS) results and it can be seen that they **match** except by:
> _in terms of topology_
  * entries HHH (Hie-Hie-Hie), JJJ (Hip-Hip-Hip), OOO (Hid-Hid-Hid), RRR (Arg-Arg-Arg) and WWW (Trp-Trp-Trp) have some improper dihedrals inverted;
  * WWW which has 3 extra improper dihedrals related to atoms sharing in the 5-ring and 6-ring of the TRP. These improper dihedrals would be there in order to keep the planarity between 5-ring and 6-ring;
  * YYY (Tyr-Tyr-Tyr) for which atom CZ (id 15, 36 and 57) got atom type CA instead of C in GROMACS.
> _in terms of parameters_
  * charges (can be either _gasteiger_ or _bcc_ with ACPYPE);
  * YYY: 6 bonds and 9 dihedrals, all involving atom CZ (id 15, 36 and 57), because of atom type changing mentioned above;

When we look the bonded potential energies for all entries, the difference is no bigger than 0.002% (even for the entries with inverted improper dihedrals) with one solely exception: YYY with 1.9%, because of the atom type change seen above.

The source code for these tests can be seen [here](http://code.google.com/p/acpype/source/browse/trunk/test/check_acpype.py).

See detailed output [here](ResiduesTest.md)

## ACPYPE and OPLS/AA Parameters Generation (Experimental) ##

The biggest problem here is to guess the correct atom type according to OPLS/AA
forcefield.

Looking at `.../gromacs/top/ffoplsaanb.itp` file and one can see e.g.:
```
[ atomtypes ]
; full atom descriptions are available in ffoplsaa.atp
; name    bond_type    mass         charge   ptype   sigma        epsilon
...
 opls_015   C2  6      14.02700     0.285       A    3.80000e-01  4.93712e-01 ; SIG
 opls_016   C2  6      14.02700    -0.100       A    3.90500e-01  4.93712e-01 ; SIG
...
```

In this example, atom types `opls_015` and `opls_016` shared the same `bond_type` but differs
by charge and van der waals parameters. Having the same `bond_type` means that
they share the same _bonded_ parameters as seen in `.../gromacs/top/ffoplsaabon.itp`.

Since ACPYPE uses ANTECHAMBER, for it, all atom types found will be also the `bond_type`
and their charges will be calculated specifically. So the "trick" here is to map the
GAFF or AMBER atom type against OPLS atom type.

I did it by creating a mapping dictionary that links GAFF or AMBER atom type to a list
of possible OPLS atom types based on the 22 pdb entries mentioned above by comparing
GROMACS' OPLS output with ACPYPE's GAFF and AMBER results.

This is not ideal. The correct approach would involve determining the atom chemical function
based on its surrounding chemical neighbourhood, which ANTECHAMBER does but not for
OPLS/AA atom types.

Nevertheless, I hope ACPYPE can be helpful since it also shows in the ACPYPE
GROMACS OPLS _itp_ output file a suggestion of possible replacements for the OPLS
atom types guessed by ACPYPE.

Another alternative would be to port the algorithm developed in
[MKTOP](http://labmm.iq.ufrj.br/mktop) (done in _perl_ and following the correct
approach mentioned above) that determines the OPLS/AA atom types to _python_ and
then add it to ACPYPE (under study).

Still, since that for `bond_types`, the chances that one picked the right `bond_type`
are almost sure (but not always, sometimes a CT can be CT\_2 or CT\_3 in OPLS/AA),
the _bonded_ parameters (i.e., the topology) determined by ACPYPE are fundamentally
correct, leaving just one opened end: the _bonded_ parameters may not be present
in `.../gromacs/top/ffoplsaabon.itp`, and for that `grompp` will fail showing which
lines in your _itp_ file have missing parameters.

However, calculated parameters derived for AMBER99SB/GAFF are there and one can just
enable them by commenting out (remove "`;`") in that line.

**BUT**, although many parameters for AMBER99SB may be identical to the ones for
OPLS/AA, e.g.:
```
(from .../gromacs/top/ffamber99sbbon.itp)
X   CA  CA  X     3    30.33400     0.00000   -30.33400     0.00000     0.00000     0.00000 ; intrpol.bsd.on C6H6

(from .../gromacs/top/ffoplsaabon.itp)
  X      CA     CA     X       3     30.33400   0.00000 -30.33400   0.00000   0.00000   0.00000 ; aromatic ring
```

these parameters were NOT calculated for OPLS/AA and so one will find several divergences, e.g.:
```
(from .../gromacs/top/ffamber99sbbon.itp)
CT  CT  C   N     3     2.51040     4.18400     1.67360    -6.69440     0.00000     0.00000     ; new99sb

(from .../gromacs/top/ffoplsaabon.itp)
  CT     CT     C      N       3      4.83252  -7.65254   1.68196   1.13805   0.00000   0.00000 ; propanamide
```

As for improper dihedrals, one can assume they are if not identical, at least equivalent, so one can
use the parameters derived for AMBER99SB by ACPYPE/ANTECHAMBER with OPLS/AA
by removing "`;`" in the respective lines of your generated _itp_ file.

Nevertheless, an ideal solution would involve finding the right OPLS atom types and parameters
calculation according to [Jorgensen et al. (1996)](http://dx.doi.org/10.1021%2Fja9621760).
Be aware of it!

By the way, it is important to emphasise that the mechanism used by ANTECHAMBER to
find the topology (i.e., how atoms are connected in bonds, angles and dihedrals)
of a chemical compound is pretty reliable for an all-atoms forcefields-like (it is just
the way how AMBER forcefields are developed).

And a last observation. If one finds in the GROMACS OPLS _itp_ file generated an
atom type like `opls_x` with mass `0.000`, it's because my mapping dictionary failed
to guess a putative OPLS atom type. Please let me know about it if you run into it.

## References ##

  * **[Antechamber](http://ambermd.org/#AmberTools)** ([manual](http://ambermd.org/doc11/AmberTools.pdf))
    1. J. Wang, W. Wang, P.A. Kollman and D.A. Case. "Automatic atom type and bond type perception in molecular mechanical calculations". Journal of Molecular Graphics and Modelling, 25, 247-260 (2006).
    1. [J. Wang, R.M. Wolf, J.W. Caldwell, P.A. Kollman and D.A. Case. "Development and testing of a general AMBER force field". Journal of Computational Chemistry, 25, 1157-1174 (2004)](http://ambermd.org/antechamber/gaff.pdf).
  * **[AMBER](http://ambermd.org/#Amber11)** ([manual](http://ambermd.org/doc11/Amber11.pdf))
    1. [D.A. Case, T.E. Cheatham, III, T. Darden, H. Gohlke, R. Luo, K.M. Merz, Jr., A. Onufriev, C. Simmerling, B. Wang and R. Woods. "The Amber biomolecular simulation programs". J. Computat. Chem. 26, 1668-1688 (2005)](http://people.cs.vt.edu/~onufriev/PUBLICATIONS/amber8JCC.pdf).
    1. [J.W. Ponder and D.A. Case. "Force fields for protein simulations". Adv. Prot. Chem. 66, 27-85 (2003)](http://dasher.wustl.edu/bio5476/reading/advprotchem-66-27-03.pdf).
  * **[GROMACS](http://www.gromacs.org)** ([manual](http://www.gromacs.org/@api/deki/files/82/=gromacs4_manual.pdf))
    1. B. Hess, C. Kutzner, D. van der Spoel and E. Lindahl. "GROMACS 4: Algorithms for Highly Efficient, Load-Balanced, and Scalable Molecular Simulation". J. Chem. Theory Comput., 4, 3, 435-447 (2008).
  * **[ffAMBER](http://ffamber.cnsm.csulb.edu)**
    1. [Eric J. Sorin and Vijay S. Pande. "Exploring the Helix-Coil Transition via All-Atom Equilibrium Ensemble Simulations". Biophysical Journal, 88 (4), 2472-2493 (2005)](http://ffamber.cnsm.csulb.edu/pdfs/sorin_amber99phi_2005bj.pdf).
    * AMBER-94 => ffamber94: [Cornell et al. (1995), JACS 117, 5179-5197](http://ffamber.cnsm.csulb.edu/pdfs/cornell_amber94_1995jacs.pdf)
    * AMBER-96 => ffamber96: [Kollman (1996), Acc. Chem. Res. 29, 461-469](http://ffamber.cnsm.csulb.edu/pdfs/kollman_amber_1996acr.pdf)
    * AMBER-99 => ffamber99: [Wang et al. (2000), J. Comp. Chem. 21, 1049-1074](http://ffamber.cnsm.csulb.edu/pdfs/wang_amber99_2000jcc.pdf)
    * AMBER-03 => ffamber03: [Duan et al. (2003), J. Comp. Chem. 24, 1999-2012](http://ffamber.cnsm.csulb.edu/pdfs/duan_AMBER03_2003jcc.pdf)
    * AMBER-GS => ffamberGS: [Garcia & Sanbonmatsu (2002), PNAS 99, 2782-2787](http://ffamber.cnsm.csulb.edu/pdfs/garcia_amberGS_2002pnas.pdf)
    * AMBER-GS-S => ffamberGSs: [Nymeyer & Garcia (2003) PNAS 100, 13934-13939](http://ffamber.cnsm.csulb.edu/pdfs/garcia_amberGS-S_2003pnas.pdf)
    * AMBER-99f => ffamber99p: [Sorin & Pande (2005). Biophys. J. 88(4), 2472-2493](http://ffamber.cnsm.csulb.edu/pdfs/sorin_amber99phi_2005bj.pdf)
    * AMBER-99SB => ffamber99sb: [Hornak et. al (2006). Proteins 65, 712-725](http://ffamber.cnsm.csulb.edu/pdfs/simmerling_amber99sb_2006proteins.pdf)
    * AMBER-99BSC0 (not ported yet): Pérez et al. (2007). Biophys J 92, 3817–29
  * **[CNS](http://cns.csb.yale.edu)**
    1. [A.T. Brunger, P.D. Adams, G.M. Clore, P.Gros, R.W. Grosse-Kunstleve, J.-S. Jiang, J. Kuszewski, N. Nilges, N.S. Pannu, R.J. Read, L.M. Rice, T. Simonson, G.L. Warren. "Crystallography & NMR System (CNS), A new software suite for macromolecular structure determination". Acta Cryst. D54, 905-921 (1998)](http://cns.csb.yale.edu/v1.2/about_cns/jn0043.pdf).
    1. [A.T. Brunger. "Version 1.2 of the Crystallography and NMR System". Nature Protocols 2, 2728-2733 (2007)](http://cns.csb.yale.edu/v1.2/about_cns/brunger_nature_protocols_2007.pdf).
  * **[XPLOR-HIH](http://nmr.cit.nih.gov/xplor-nih/)**
    1. C.D. Schwieters, J.J. Kuszewski, N. Tjandra and G.M. Clore. "The Xplor-NIH NMR Molecular Structure Determination Package". J. Magn. Res., 160, 66-74 (2003).
    1. C.D. Schwieters, J.J. Kuszewski and G.M. Clore. "Using Xplor-NIH for NMR molecular structure determination". Progr. NMR Spectroscopy 48, 47-62 (2006).
  * **CNS/XPLOR-NIH:** For a discussion of the effects of using the Amber/GB force field (rather than the Engh-Huber-based targets) see B. Xia, V. Tsui, D.A. Case, H.J. Dyson, and P.E. Wright. "Comparison of protein solution structures refined by molecular dynamics simulation in vacuum, with a generalized Born model, and with explicit water". J. Biomol. NMR 22, 317-331 (2002).
  * **[XPLO2D](http://xray.bmc.uu.se/usf/xplo2d_man.html)** ([FAQ](http://xray.bmc.uu.se/hicup/faq.html), [download](ftp://xray.bmc.uu.se/pub/gerard/xutil/))
    1. G.J. Kleywegt, K. Henrick, E.J. Dodson and D.M.F. van Aalten. "Pound-wise but penny-foolish: How well do micromolecules fare in macromolecular refinement?" Structure, 11, 1051-1059 (2003).
  * **[CHARMM](http://www.charmm.org)**
    1. B.R. Brooks, R.E. Bruccoleri, B.D. Olafson, D.J. States, S. Swaminathan and M. Karplus. "CHARMM: A Program for Macromolecular Energy, Minimization, and Dynamics Calculations". J. Comp. Chem. 4, 187-217 (1983).
  * **[NAMD](http://www.ks.uiuc.edu/Research/namd)** ([manual](http://www.ks.uiuc.edu/Research/namd/current/ug.pdf))
    1. [James C. Phillips, Rosemary Braun, Wei Wang, James Gumbart, Emad Tajkhorshid, Elizabeth Villa, Christophe Chipot, Robert D. Skeel, Laxmikant Kale and Klaus Schulten. "Scalable molecular dynamics with NAMD". Journal of Computational Chemistry. 26, 1781-180 (2005)](http://www.ks.uiuc.edu/Publications/Papers/PDF/PHIL2005/PHIL2005.pdf).
  * **[VMD](http://www.ks.uiuc.edu/Research/vmd)** ([manual](http://www.ks.uiuc.edu/Research/vmd/current/ug.pdf))
    1. [W. Humphrey, A. Dalke and K. Schulten. "VMD - Visual Molecular Dynamics". J. Molec. Graphics. 14, 33-38 (1996)](http://www.ks.uiuc.edu/Publications/Papers/PDF/HUMP96/HUMP96.pdf).
  * **[Empirical Force Field Parameters](http://www.psc.edu/general/software/packages/charmm/tutorial/mackerell/PARAM_00.pdf)** by Alexander MacKerell
  * **[YASARA Autosmiles](http://www.yasara.org/autosmiles.htm)**
  * **[OPLS](http://en.wikipedia.org/wiki/OPLS)**
    1. [W.L. Jorgensen, D.S. Maxwell and J. Tirado-Rives. "Development and Testing of the OPLS All-Atom Force Field on Conformational Energetics and Properties of Organic Liquids". J. Am. Chem. Soc. 118: 11225–11236 (1996)](http://dx.doi.org/10.1021%2Fja9621760).