import argparse
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
import random
import pandas as pd
from collections import Counter
random.seed(10)

"""
This script simulates the decomposition of Milled Wood Lignin (MWL) into Deep Eutectic Solvent Lignin (DESL) fragments 
and calculates various molecular characteristics for both the MWL and DESL structures. 

Key functionalities:
1. **SMILES Representation**: Configures MWL structure using its SMILES representation.
2. **Molecular Characteristics**: Computes molecular weight, functional group counts, and poly dispersity index (PDI) for both the MWL and DESL structures.
3. **DESL Simulation**: Simulates the decomposition of MWL into DESL fragments using the `simulate_DES_deployermization` function.

Imports:
- `rdkit`: Used for manipulating molecular structures (SMILES, molecular weight, etc.).
- `utils.desl_simulation`: Contains custom functions for finding functional group counts, calculating PDI, and simulating the decomposition of MWL.

This script is designed to help in the analysis and visualization of the molecular transformations involved in DESL formation from MWL.

Sample command:
python utils/desl_simulation.py --ether-cleavage 1.0 --image-out image/desl_output.png

Parameter usage:
- --smiles: Input MWL structure in SMILES format. If omitted, notebook default is used.
- --ether-cleavage: Fraction of ether bonds to cleave (0.0 to 1.0).
- --image-out: Optional image file path to save MWL vs DESL visualization.

Author: Sudha Eswaran
Date: 03/15/2024
"""


DEFAULT_MWL_SMILES = "OC=1C=CC(=CC1OC)C(O)C(OC=2C=CC(=CC2OC)C(O)C(OC3=C(OC)C=C(C=C3C4=CC(=CC(OC)=C4OC(CO)C(O)C=5C=CC(OC(CO)C(O)C=6C=CC(O)=C(OC)C6)=C(OC)C5)C7OCC8C(OCC78)C9=CC(OC)=C%10OC(C%11=CC(OC)=C%12OC(C=%13C=CC(O)=C(OC)C%13)C(C%12=C%11)CO)C(C%10=C9)CO)C(O)C(OC=%14C=CC(=CC%14OC)C(O)C(OC=%15C=CC(=CC%15OC)C(O)C(OC=%16C=CC(=CC%16OC)C%17OCC%18C(OCC%17%18)C=%19C=CC(O)=C(OC)C%19)CO)CO)CO)CO)CO"


def simulate_DES_deployermization(smile, ethercleavage):
    """
    Simulates the decomposition of a molecule by cleaving C-O-C bonds and generating ketones.

    Args:
        smile (str): The SMILES string of the molecule to simulate.
        ethercleavage (float): The percentage of C-O-C bonds to cleave.

    Returns:
        str: The SMILES representation of the decomposed and modified molecule.
    """
    mol = Chem.MolFromSmiles(smile)
    mols = fragment(mol, ethercleavage)
    return generate_ketones(mols)


def find_functional_groups_counts(molsmile, *fntype):
    """
    Identifies functional groups in a molecule based on its SMILES representation and counts their occurrences.

    Args:
        molsmile (str): The SMILES string of the molecule to analyze.
        *fntype (str): Optional. Specific functional group type to count. If not provided, counts all functional groups.

    Returns:
        dict or int: If no functional group type is specified, returns a dictionary with the counts of various functional groups.
                      If a functional group type is specified, returns the count of that specific functional group.
    """
    fntype = fntype[0] if len(fntype) != 0 else ""
    functionalgroups = {'ether': '[C;H1][OX2][c]',
                        'carbonyl': '[CX3]=[OX1]',
                        'methoxy': '[c][OX2][CH3]',
                        'phenolicOH': 'c1ccccc1[OH]',
                        'aliphaticOH_primary': 'cCCC[OH]',
                        'aliphaticOH_secondary': 'cC[OH]',
                        'cinnamyl_alcohol': 'c1ccccc1C=CCO'}
    mol = Chem.MolFromSmiles(molsmile)
    ether = mol.GetSubstructMatches(
        Chem.MolFromSmarts(functionalgroups['ether']))
    carbonyl = mol.GetSubstructMatches(
        Chem.MolFromSmarts(functionalgroups['carbonyl']))
    methoxy = mol.GetSubstructMatches(
        Chem.MolFromSmarts(functionalgroups['methoxy']))
    phenol = mol.GetSubstructMatches(
        Chem.MolFromSmarts(functionalgroups['phenolicOH']))
    aliphaticOH_secondary = mol.GetSubstructMatches(
        Chem.MolFromSmarts(functionalgroups['aliphaticOH_secondary']))
    aliphaticOH_primary = mol.GetSubstructMatches(
        Chem.MolFromSmarts(functionalgroups['aliphaticOH_primary']))
    cinnamyl_alcohol = mol.GetSubstructMatches(
        Chem.MolFromSmarts(functionalgroups['cinnamyl_alcohol']))

    funcgrp_dict = {
        "ether": len(ether),
        "carbonyl": len(carbonyl),
        "methoxy": len(methoxy),
        "phenolicOH": len(phenol),
        "aliphaticOH_secondary": len(aliphaticOH_secondary),
        "aliphaticOH_primary": len(aliphaticOH_primary),
        "cinnamyl_alcohol": len(cinnamyl_alcohol),
    }

    if fntype in funcgrp_dict.keys():
        return funcgrp_dict[fntype]
    else:
        return funcgrp_dict


def mol_with_atom_index(mol):
    """
    Adds atom indices to a molecule for visualization purposes.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object.

    Returns:
        rdkit.Chem.Mol: The input molecule with atom indices added.
    """
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol


def fragment_on_bond(mol, atom1, atom2):
    """
    Fragments a molecule at the specified bond between two atoms.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object to be fragmented.
        atom1 (int): The index of the first atom in the bond.
        atom2 (int): The index of the second atom in the bond.

    Returns:
        rdkit.Chem.Mol: The fragmented molecule.
    """
    bond = mol.GetBondBetweenAtoms(atom1, atom2)
    new_mol = mol
    if bond is not None:
        new_mol = Chem.FragmentOnBonds(
            mol, [bond.GetIdx()], dummyLabels=[(4, 4)])
    return new_mol


def fragment(mol, cleave_per):
    """
    Fragments a molecule by cleaving C-O-C bonds based on a given cleavage percentage.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object to be fragmented.
        cleave_per (float): The percentage of bonds to cleave, between 0 and 1.

    Returns:
        rdkit.Chem.Mol: The fragmented molecule.
    """
    ether_bond_tuples = mol.GetSubstructMatches(
        Chem.MolFromSmarts("[C;H1,H2][OX2][c]"))
    cnt = int(len(ether_bond_tuples) * cleave_per)
    new_mol = mol

    for frgs in random.sample(ether_bond_tuples, k=len(ether_bond_tuples)):
        i = 2
        while i < len(frgs) and cnt > 0:
            a = mol.GetAtomWithIdx(frgs[i-2]).GetAtomicNum()
            b = mol.GetAtomWithIdx(frgs[i-1]).GetAtomicNum()
            c = mol.GetAtomWithIdx(frgs[i]).GetAtomicNum()
            if a == 6 and b == 8 and c == 6:
                new_mol = fragment_on_bond(new_mol, frgs[i-2], frgs[i-1])
            i += 1
            cnt -= 1
    return update_phenolic_hydroxyl(new_mol)


def update_phenolic_hydroxyl(mol):
    """
    Updates phenolic hydroxyl groups by replacing them with hydroxyl groups.

    Args:
        mol (rdkit.Chem.Mol): The RDKit molecule object to be modified.

    Returns:
        rdkit.Chem.Mol: The modified molecule with updated phenolic hydroxyl groups.
    """
    mol.GetSubstructMatches(Chem.MolFromSmarts("[4*]O"))
    mod_mol = Chem.ReplaceSubstructs(mol, Chem.MolFromSmiles(
        'O[4*]'), Chem.MolFromSmiles('O[H]'), replaceAll=True)
    return mod_mol[0]


def generate_ketones(mols):
    """
    Generates ketones from a list of molecules by modifying specific structures.

    Args:
        mols (rdkit.Chem.Mol): The RDKit molecule(s) to generate ketones from.

    Returns:
        str: The SMILES representation of the generated ketones.
    """
    ketone_type_list = ['vanillin', 'HB1', 'HB2']
    mol_smile = Chem.MolToSmiles(mols)
    mol_list = mol_smile.split(".")
    upd_mol_list = []
    for cmp in mol_list:
        mol = Chem.MolFromSmiles(cmp)
        act_site = mol.GetSubstructMatches(Chem.MolFromSmarts('[4*][C][c]'))
        if len(act_site) > 0:
            em0 = Chem.EditableMol(mol)
            em0.RemoveAtom(act_site[0][0])
            mol = em0.GetMol()
        ss = mol.GetSubstructMatches(Chem.MolFromSmarts('[OH][C][c]'))
        if len(ss) > 0:
            mod_mol = Chem.ReplaceSubstructs(mol, Chem.MolFromSmiles(
                'C[4*]'), Chem.MolFromSmiles('C=O'), replaceAll=True)
            mod_mol2 = mod_mol[0]
        else:
            mod_mol = Chem.ReplaceSubstructs(mol, Chem.MolFromSmiles(
                'C[4*]'), Chem.MolFromSmiles('C'), replaceAll=True)
            mod_mol2 = mod_mol[0]
        ss = mod_mol2.GetSubstructMatches(Chem.MolFromSmarts('[OH][C][c]'))
        while len(ss) > 0:
            em1 = Chem.EditableMol(mod_mol2)
            em1.RemoveAtom(ss[0][0])
            mod_mol2 = em1.GetMol()
            ss = mod_mol2.GetSubstructMatches(Chem.MolFromSmarts('[OH][C][c]'))

        ketone_type = random.choice(ketone_type_list)
        if ketone_type == 'HB2':
            ss2 = mod_mol2.GetSubstructMatches(
                Chem.MolFromSmarts('[OH][C][C]'))
            if len(ss2) > 0:
                em1 = Chem.EditableMol(mod_mol2)
                em1.RemoveAtom(ss2[0][0])
                mod_mol2 = em1.GetMol()
        upd_mol_list.append(Chem.MolToSmiles(mod_mol2))
    return ".".join(upd_mol_list)


def calculate_pdi(smile):
    """
    Calculates the Polydispersity Index (PDI) of a molecule based on its SMILES string.

    Args:
        smile (str): The SMILES string of the molecule.

    Returns:
        tuple: The weight-average molecular weight, number-average molecular weight, and PDI.
    """
    compounds = smile.split(".")
    mol_wt_compound = Descriptors.ExactMolWt(Chem.MolFromSmiles(smile))
    nA_Mw_sum = 0
    wt_Mw_sum = 0
    canon_smile = [Chem.CanonSmiles(x) for x in compounds]
    canon_smile_count = Counter(canon_smile)
    dict_cnt = dict(canon_smile_count)
    for molsmile, cnt in dict_cnt.items():
        mol_wt_frag = Descriptors.ExactMolWt(Chem.MolFromSmiles(molsmile))
        nA_Mw_sum += mol_wt_frag * cnt
    for molsmile, cnt in dict_cnt.items():
        mol_wt_frag = Descriptors.ExactMolWt(Chem.MolFromSmiles(molsmile))
        mol_wt_fraction = (mol_wt_frag * cnt) / nA_Mw_sum
        wt_Mw_sum += mol_wt_frag * mol_wt_fraction
    nA_Mw = nA_Mw_sum / sum(dict_cnt.values())
    return wt_Mw_sum, nA_Mw, wt_Mw_sum / nA_Mw


def calculate_dp(smile):
    """
    Calculates the degree of polymerization (DP) from the SMILES string of a molecule.

    Args:
        smile (str): The SMILES string of the molecule.

    Returns:
        list: A list of unique degree of polymerization values for the molecule.
    """
    compounds = smile.split(".")
    dp = []
    for molsmile in compounds:
        dp.append(int(Descriptors.ExactMolWt(
            Chem.MolFromSmiles(molsmile)) / 180))
    return list(set(dp))


def calculate_molecular_weight(smiles):
    """
    Calculates the exact molecular weight of a molecule based on its SMILES string.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        float or None: The exact molecular weight if the SMILES is valid, otherwise None.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return Descriptors.ExactMolWt(mol)
    else:
        return None


def build_parser():
    """Builds a command-line parser for running DESL simulation as a script."""
    parser = argparse.ArgumentParser(
        description=(
            "Simulate decomposition of MWL into DESL fragments and print molecular metrics."
        )
    )
    # Input polymer structure in SMILES form.
    parser.add_argument(
        "--smiles",
        default=DEFAULT_MWL_SMILES,
        help="Input MWL SMILES string (defaults to the notebook example).",
    )
    # Controls how aggressively ether linkages are broken during simulation.
    parser.add_argument(
        "--ether-cleavage",
        type=float,
        default=1.0,
        help="Fraction of ether bonds to cleave (0.0 to 1.0). Default: 1.0",
    )
    # Optional output image path for side-by-side molecular visualization.
    parser.add_argument(
        "--image-out",
        default="",
        help="Optional path to save a comparison image (MWL vs DESL).",
    )
    return parser


def main():
    """Runs the same simulation workflow used in desl-simulation.ipynb."""
    parser = build_parser()
    args = parser.parse_args()

    if not 0.0 <= args.ether_cleavage <= 1.0:
        raise ValueError("--ether-cleavage must be between 0.0 and 1.0")

    smile = args.smiles
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        raise ValueError("Invalid input SMILES provided via --smiles")

    canonical_smiles = Chem.CanonSmiles(smile)
    molecular_weight = int(Descriptors.ExactMolWt(mol))
    functional_groups = find_functional_groups_counts(smile)
    pdi = round(calculate_pdi(smile)[2], 1)

    print("\n\nDES Lignin : Simulation Model")
    print("=============================\n")

    print("\nINPUT")
    print(f"Molecular structure of MWL (SMILES) : \n{canonical_smiles}")
    print(f"\nMolecular weight (g/mol)           : {molecular_weight}")
    print(f"\nFunctional group                   : {functional_groups}")
    print(f"\nPoly dispersity index (PDI)        : {pdi}")

    smi = simulate_DES_deployermization(smile, args.ether_cleavage)
    mol_smi = Chem.MolFromSmiles(smi)

    print("\nOUTPUT")
    print(f"\nNumber of fragments  : {len(smi.split('.'))}")
    print("\t" + "\n\t".join(smi.split('.')))
    print(
        f"\nTotal Molecular weight (g/mol) of all fragments: {int(Descriptors.ExactMolWt(mol_smi))}"
    )
    print(
        f"\nFunctional group                                 : {find_functional_groups_counts(smi)}"
    )
    print(
        f"\nPoly dispersity index (PDI)                      : {round(calculate_pdi(smi)[2], 1)}"
    )

    if args.image_out:
        labels = [
            "Molecular structure of MWL",
            "Molecular structure after simulation (DESL)",
        ]
        grid = Draw.MolsToGridImage(
            [mol, mol_smi],
            molsPerRow=2,
            subImgSize=(600, 600),
            legends=labels,
        )
        grid.save(args.image_out)
        print(f"\nSaved visualization image to: {args.image_out}")


if __name__ == "__main__":
    main()
