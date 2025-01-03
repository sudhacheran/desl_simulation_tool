# DES Lignin Simulation Model

This repository contains a script that simulates the decomposition of Milled Wood Lignin (MWL) into Deep Eutectic Solvent Lignin (DESL) fragments. The script also computes various molecular characteristics for both the MWL and DESL structures, including molecular weight, functional groups, and poly dispersity index (PDI). Additionally, it provides molecular visualizations to aid in understanding the decomposition process.
Features

**Molecular Simulation:** Decomposes MWL structure into DESL fragments using a simulation function.
**Molecular Characteristics:** Computes important molecular properties like molecular weight, functional groups, and poly dispersity index (PDI).
**SMILES Representation:** Takes MWL as input in SMILES format and outputs the decomposed DESL structure in SMILES format.
**Modular Functions:** Custom functions for functional group counting, PDI calculation, and simulation of decomposition.
**Visualization:** Draws molecular structures of MWL and DESL using RDKit.

## **Requirements**
Python 3.x
Required libraries:
rdkit: For molecular manipulation and visualization.
utils.desl_simulation: A custom module containing functions for simulation and analysis.

To install the required libraries, you can use pip:
pip install -r requirements.txt

## **Script Overview:**
The script begins by configuring the MWL structure in SMILES format.
It then computes important molecular characteristics such as molecular weight, functional group counts, and poly dispersity index (PDI) for the MWL structure.
The script simulates the decomposition of MWL into DESL fragments and calculates properties of the fragments.
Finally, it visualizes the MWL and DESL structures using RDKit's drawing functions.

**Running the script:**
You can use the attached desl-simulation.ipynb Jupyter Notebook to run the script interactively. The notebook allows you to execute the script step-by-step in a more interactive environment.
Alternatively, you can run the script directly in a Python environment by executing:
`python your_script.py`
This will display the molecular characteristics of both MWL and DESL, as well as the visualization of the molecular structures.

Input:
The input for the script is the SMILES string representing the MWL structure. You can modify the smile variable in the script with any valid SMILES string to simulate different lignin structures.
Output:
The script outputs the following:
SMILES string of the MWL and DESL structures.
Molecular weight of both MWL and DESL.
Count of functional groups in MWL and DESL.
Poly dispersity index (PDI) for MWL and DESL.
Visualization of MWL and DESL structures displayed as a grid image.

**Example Output**

Input:
SMILES string as input for the MWL structure:
OC=1C=CC(=CC1OC)C(O)C(OC=2C=CC(=CC2OC)C(O)C(OC=3C=CC(=CC3OC)C(O)C(OC=4C=CC(=CC4OC)C(O)C(OC=5C=CC(=CC5OC)C(O)C(OC=6C=CC(=CC6OC)C7OCC8C(OCC78)C9=CC(OC)=C%10OC(C%11=CC(OC)=C%12OC(C=%13C=CC(O)=C(OC)C%13)C(C%12=C%11)CO)C(OC%14=C(OC)C=C(C=C%14C%10=C9)C%15OCC%16C(OCC%15%16)C%17=CC(OC)=C%18OC(C=%19C=CC(OC(CO)C(O)C=%20C=CC(O)=C(OC)C%20)=C(OC)C%19)C(C%18=C%17)CO)CO)CO)CO)CO)CO)CO

Output:
After running the script, the output will display information such as:
DES Lignin : Simulation Model
=============================

INPUT
======
Molecular structure of MWL (SMILES) :
OC=1C=CC(=CC1OC)C(O)C(OC=2C=CC(=CC2OC)C(O)C(OC=3C=CC(=CC3OC)C(O)C(OC=4C=CC(=CC4OC)C(O)C(OC=5C=CC(=CC5OC)C(O)C(OC=6C=CC(=CC6OC)C7OCC8C(OCC78)C9=CC(OC)=C%10OC(C%11=CC(OC)=C%12OC(C=%13C=CC(O)=C(OC)C%13)C(C%12=C%11)CO)C(OC%14=C(OC)C=C(C=C%14C%10=C9)C%15OCC%16C(OCC%15%16)C%17=CC(OC)=C%18OC(C=%19C=CC(OC(CO)C(O)C=%20C=CC(O)=C(OC)C%20)=C(OC)C%19)C(C%18=C%17)CO)CO)CO)CO)CO)CO)CO

Molecular weight (g/mol)           : 3520
Functional group                   : {'phenolic': 6, 'ether': 10}
Poly dispersity index (PDI)        : 2.1

OUTPUT
======

Number of fragments (DESL Structures)  : 3
	OC=1C=CC(=CC1OC)C(O)
	OC=2C=CC(=CC2OC)C(O)
	OC=3C=CC(=CC3OC)C(O)

Total Molecular weight (g/mol) of all fragments : 980
Functional group                                 : {'phenolic': 2, 'ether': 4}
Poly dispersity index (PDI)                      : 1.5

**Visualization:**
The script also generates and displays a grid image of the MWL and DESL molecular structures.
Contributing

<hr/>
Feel free to fork this repository, submit issues, or send pull requests to improve the simulation and the overall functionality of the script.

License

