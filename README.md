Script: `alchemical_setup.py`

An automatic setup of relative free energy calculations.

The tool generates input files needed to perform relative solvation
and binding free energy calculations with GROMACS.

To provide more details on the algorithmic side, we achieve the above by performing the following steps:

* Initialize two objects, one for each molecule, A and B (apart from containing information on the `.top` and `.gro` file names, one of its attributes will be reserved for a dictionary to store the atom indices mapping)
* Initialize an object for chimeric molecule whose attributes are
	* a list keep track of the types of dummy atoms,
	* a dictionary to contain atom types of the atoms with new indices,
	* a nested dictionary with its lists to store entries for the `[bonds]`, `[angles]`, and `[dihedrals]` directives,
	* a nested dictionary with its dictionaries to store atom type sequences)
* Identify the atom indices of each molecule in MCSS by reading in the `map.txt` file
* Build the `[atoms]` directive read in and store appropriately all the fields of the `[atoms]` directive for both `.top` files to avoid duplicates, skip the MCSS atoms of molecule B assign new atom names sort atoms by their new name
* Build the `[atomtypes]` directive
* Build the `[pairs]` directive
* Build directives for the bonded interactions
* Prepare the final `.gro` file
	* extract the coordinates of atoms of both molecules
	* find the coordinates of atoms of molecule B when its MCSS overlaid on that of molecule A
	* append coordinates of the peripheral atoms of molecule B to the molecule A `.gro` file
	* renumber entries in the `.gro` file
* Write out the final `.top` file



Help for `alchemical_setup.py` (obtained with `python alchemical_setup.py -h`) is:

```Options:
  -h, --help  show this help message and exit
  -m MAPTXT   The atom map .txt file.
  -a TOP_A    The .top file of molecule A.
  -b TOP_B    The .top file of molecule B.
  -s MCSS     The substructure .mol2 file.
  -o OUT_DIR  The output directory name.
  -v VERBOSE  Verbosity.
```
