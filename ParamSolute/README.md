# Parametrize small molecules

In this tutorial I will explain how _I_ parametrize small molecules, e.g., when I have series of reasonably small molecules that I want to simulate. I will use the general Amber force field (GAFF) together with AM1-BCC partial charges for the solutes.

You will  need a range of in-house Python scripts that I have created and that you find in this Github repository: [Scripts](http://www.github.com/sgenheden/scripts). Download them to your local machine. I will be referring to this folder with `$SCRIPTS`.

Also make sure that you have the prerequisites python libraries.

Finally, you also need to install the [Amber tools](http://www.ambermd.org/) and make sure that the executables are in the system path.

 If there is just one solute, the procedure goes something like this
1. Create / find a structure of the solute
2. Run the structure through `antechamber` to obtain topology, GAFF atom types and AM1-BCC charges
3. Run the solue through `parmchk` to guess any missing bonded parameters
4. Put together a parameter + topology file in Amber format using `tleap`
5. Optionally, convert the **prmtop** file to either Gromacs or Lammps format

and if you only have one or a few solutes, you could improve upon this recipe by computing RESP charges or thoroughly parametrize missing bonded parameters. However, for most applications the standard parameters will be alright.

However, if you have a lot of solutes, the procedure above becomes tedious. Therefore, I will show you how to use a range of in-house scripts that automatize this procedure.

## Obtaining 3D structures

In this tutorial we will parameterize benzene, phenol and toluene. The first task is to obtain the structure, which will be accomplished by
1. Obtain the _SMILES_ strings of the solutes.
2. Send the SMILES strings to a web-server.

The easiest way to find the SMILES is to go to the [ChemSpider](http://www.chemspider.com/) web-service. Type in the name of a solute and under _More details:_ you will find the SMILES.

For our solutes it should be

    benzene: c1ccccc1
    phenol: c1ccc(cc1)O
    toluene: Cc1ccccc1

Let's put them in a file together

    cat << EOF > smiles.txt
    benzene c1ccccc1
    phenol c1ccc(cc1)O
    toluene Cc1ccccc1
    EOF

Now when we have obtained the SMILES strings, we can use the following command to obtain the structure in an xyz-file.

    python $SCRIPTS/Mol/smiles2xyz.py --inlist smiles.txt


**Using the ChemSpider python interface**

As an alternative to the procedure above, we can use the excellent python interface to the ChemSpider database. You first have to register on the ChemSpider homepage so that you can have a private security token. You can read more about it [here](http://www.chemspider.com/aboutservices.aspx), but simply put: you need to register and then you can see your token on the *Profile* page.

The web services that ChemSpider provides are nicely wrapped in a python package [ChemSpiPy](http://chemspipy.readthedocs.org/). This is most easily installed with pip.

    pip install chemspipy

Your personal token is required to setup the python environment. I usually store in the environment variable `$SPIDERKEY`. If you do that, you can use the following scripts to obtain all the SMILES directly from ChemSpider without typing in solutes manually *and* obtain the 3D-structures.

**NB!** *This might fail!* and then you have to resort to the manual procedure.

First, we put the names of all the solutes in a list

    cat << EOF > solutes.txt
    benzene
    phenol
    toluene
    EOF

Second, we obtain the structures

    python $SCRIPTS/Mol/build_mols.py solutes.txt

## Parameterizing the solutes

Now when we have the structure of all the solutes, we can parametrize them using the procedure outlined above. Fortunately, we have an in-house script that performs all of those steps automatically, for all solutes.

First, we put the names of all the solutes in a list

    cat << EOF > solutes.txt
    benzene
    phenol
    toluene
    EOF

and then we run the following command

    python $SCRIPTS/Mol/param_solutes.py solutes.txt

For each solute, you have obtained a number of files. We will take benzene as an example

    benzene.crd         the coordinates in Amber format   
    benzene.frcmod      the additional parameters, guess by parmchk
    benzene.prepi       the topology, atom types and charges in Amber format
    benzene_amb.top     the Amber parameter-topology file
    benzene_gmx.top     the Gromacs parameter-topology file
    benzene_leap.pdb    the structure in PDB format

If you want to use the solutes in Gromacs it is perhaps better to convert the .top-files to .itp-files that can be included elsewhere, in e.g. a `ffnonbonded.itp` file. As it is now the .top-files contain directives like `[ system ]` and `[ atomtypes ]` that might conflict with other such directives.

To convert, type

    python $SCRIPTS/Gromacs/top2itp.py *_gmx.top

The script will also list the atom types that need to be included elsewhere to ensure that they are defined.

    c3             6   12.01000    0.000000  A      0.339967       0.45773
    ca             6   12.01000    0.000000  A      0.339967      0.359824
    ha             1    1.00800    0.000000  A      0.259964       0.06276
    hc             1    1.00800    0.000000  A      0.264953     0.0656888
    ho             1    1.00800    0.000000  A             0             0
    oh             8   16.00000    0.000000  A      0.306647      0.880314

### Convert to Lammps format

The Lammps software is shipped with a converter from Amber format. Unfortunately, it is a little bit picky with the names of the files.

So first, rename the Amber parameter-topology files

    for X in `cat solutes.txt`
    do
    mv ${X}_amb.top ${X}.top
    done

Then we can convert all of the solutes with

    python $LMPPATH/tools/amber2lmp/amber2lammps.py

where `$LMPPATH` is the path to the folder where Lammps is installed.

and finally, rename the Amber files again, to tidy up.

    for X in `cat solutes.txt`
    do
    mv ${X}.top ${X}_amb.top
    done

The names of the Lammps data files are e.g. `data.benzene`.
