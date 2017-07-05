# Solvation free energy in water and octanol using Gromacs

In this tutorial, I will show how to compute solvation free energy of toluene in water and octanol. Both toluene and the solvents will be described with the general Amber force field (GAFF).

The MD program we will use is **Gromacs version 5.1**.

You will also need a range of in-house Python scripts that I have created and that you find in this Github repository: [Scripts](http://www.github.com/sgenheden/scripts). Download them to your local machine. I will be referring to this folder with `$SCRIPTS`.  

## Solvate the toluene molecule

We will assume that toluene already has been parametrized, so we will start directly by solvating the toluene in water and octanol using pre-equilibrated boxes.

First, we need to set the solute in a box of an appropriate size. Here, I will use a cubic box of dimensions 4 x 4 x 4 nm, which is more than sufficient for this small solute

    gmx editconf -f toluene_leap.pdb -o toluene_boxed.gro -box 4 4 4

Second, we will use the `gmx solvate` command to add the solvent molecules to this box.

For water

    gmx solvate -cp toluene_boxed.gro -cs spc216.gro -o toluene_wat.gro

and then for octanol

    gmx solvate -cp toluene_boxed.gro -cs octanol_box.gro -o toluene_oct.gro    

**NB!** please take note how many solvent molecules are added.

The `spc216.gro` box is provided by Gromacs, and although we will use TIP3P the spc box works alright. The `octanol_box.gro` is provided in this tutorial.

## Setting up the topology

Now we need to create the topologies for our systems. The topology for toluene is provided (`toluene_gmx.top`),
but we will however make it into an .itp-file so that we can use it with both solvents.

    python $SCRIPTS/Gromacs/top2itp.py toluene_gmx.top

this scripts also prints out the atom type definitions. Copy them into a new file called `ffnonboned.ff` and make sure it looks like this

    [ defaults ]    
    1               2               yes             0.5     0.833333

    [ atomtypes ]
    c3             6   12.01000    0.000000  A      0.339967       0.45773
    ca             6   12.01000    0.000000  A      0.339967      0.359824
    ha             1    1.00800    0.000000  A      0.259964       0.06276
    hc             1    1.00800    0.000000  A      0.264953     0.0656888

Next add the following atom types, after the toluene atom types

    OW           8      16.00    0.0000  A   3.15061e-01  6.36386e-01   
    HW           1       1.008   0.0000  A   0.00000e+00  0.00000e+00

which is the atom types for the TIP3P water model.

The topology file for octanol can be downloaded from [virtualchemistry.org](http://virtualchemistry.org/), it is also provided in this tutorial (`1-octanol.itp`).

Move the atom type definitions from this file to `ffnonbonded.ff`

    h1            h1      0.0000  0.0000  A   2.47135e-01  6.56888e-02
    hc            hc      0.0000  0.0000  A   2.64953e-01  6.56888e-02
    c3            c3      0.0000  0.0000  A   3.39967e-01  4.57730e-01
    oh            oh      0.0000  0.0000  A   3.06647e-01  8.80314e-01
    ho            ho      0.0000  0.0000  A   0.00000e+00  0.00000e+00

The atom types *c3* and *hc* are already provided so you don't need to copy them.

Remove the `[ defaults ]` and `[ atomtypes ]` directives and the lines that follow them. In principle the `1-octanol.itp` file should start with the `[ moleculetype ]` directive, but the comments in the beginning you can keep.

Now we can create the system topology file. For water it should look something like this and you can call it `toluene-wat.top`

    #include "ffnonbonded.itp"
    #include "toluene_gmx.itp"
    #include "amber99.ff/tip3p.itp"

    [ system ]
    Toluene in water

    [ molecules ]
    tol 1
    SOL 2160

If more or less than 2160 water molecules were added, change that number. The TIP3P topology is loaded from the Gromacs installation.

For octanol it should look something like this and you can call it `toluene-oct.top`

    #include "ffnonbonded.itp"
    #include "toluene_gmx.itp"
    #include "1-octanol.itp"

    [ system ]
    Toluene in octanol

    [ molecules ]
    tol 1
    1-octanol 174

If more or less than 174 octanol molecules were added, change that number.

## Free energy simulations

From this point on, the procedure is more or less identical for water and octanol. Therefore, I will only show the commands for octanol.

We will start with a short minimization, e.g.

    gmx grompp -f em.mdp -c toluene_wat.gro -p toluene_wat.top -o toluene_wat_em.tpr
    gmx mdrun -deffnm toluene_wat_em

For your own research project, you might want at this stage carry out some additional equilibration stage. For the purpose of this tutorial, we will skip this step.

To setup the free energy simulation, we will use the provided `dg_template.mdp` file. It will perform a 6 ns simulation per lambda state.  However, as the names implies, it is only a template and we will have to finalize it.

We will take a look at some of the free energy settings:

    free-energy              = yes
    init-lambda-state        = XXX
    coul-lambdas             = 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0
    vdw-lambdas              = 0.0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 1.0
    couple-moltype           = MMM
    couple-lambda0           = vdw-q
    couple-lambda1           = none

The first line turns on the free energy-code. The _init-lambda-state_ is the index into the lambdas arrays below. It is now set to XXX, but we will change it to 0, 1... 20. That is, we will setup 21 simulations, with different coupling parameters 0.0, 0.05 to 1.0.

The _couple-moltype_ sets which molecule type to couple with the coupling parameters, which is of course the toluene in our example. Therefore, we will change MMM too.

_couple-lambda0_ and _couple-lambda1_ instructs Gromacs to have a fully coupled solute at lambda = 0 and a fully decoupled solute at lambda = 1. For a full explanation of the settings, please look in the Gromacs manual.

To create the .tpr-files, use the following loop

    for X in {0..20}
    do
    sed -e "s/XXX/$X/" -e "s/MMM/tol/" dg_template.mdp > temp.mdp
    gmx grompp -f temp.mdp -c toluene_wat_em.gro -p toluene_wat.top -o toluene_wat_dg${X}.tpr -maxwarn 2
    done

and then run the 21 free energy simulations on a cluster.

To calculate the free energy using BAR (Bennet acceptance ratio) use the following command

    gmx bar -f toluene_wat_dg{0..20}.xvg

Remember that we compute the decoupling free energy, so the solvation free energy will be the negative of this. 
