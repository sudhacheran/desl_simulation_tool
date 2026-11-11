"""Extract the unique products from a generated DESL dataset CSV.

Each row of the dataset has a `smiles` column holding the product mixture as
dot-separated fragment SMILES. This script splits every mixture, canonicalizes
each fragment, and writes one row per distinct product with its formula, exact
mass, total number of occurrences, and number of dataset rows that contain it.

Sample command:
python data_analysis/extract_unique_products.py \
    --in data_analysis/df_data_desl_pathways_dp15_25.csv \
    --out data_analysis/unique_products_dp15_25.csv
"""

import argparse
from collections import Counter
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--in", dest="infile", required=True, help="Dataset CSV with a 'smiles' column.")
    parser.add_argument("--out", required=True, help="Output CSV of unique products.")
    args = parser.parse_args()

    df = pd.read_csv(args.infile)
    occurrences = Counter()
    rows_containing = Counter()
    invalid = 0
    for mixture in df["smiles"]:
        seen_in_row = set()
        for fragment in mixture.split("."):
            molecule = Chem.MolFromSmiles(fragment)
            if molecule is None:
                invalid += 1
                continue
            canonical = Chem.MolToSmiles(molecule)
            occurrences[canonical] += 1
            seen_in_row.add(canonical)
        rows_containing.update(seen_in_row)

    records = []
    for smiles, count in occurrences.items():
        molecule = Chem.MolFromSmiles(smiles)
        records.append(
            {
                "product_smiles": smiles,
                "formula": rdMolDescriptors.CalcMolFormula(molecule),
                "exact_mass": round(Descriptors.ExactMolWt(molecule), 3),
                "n_occurrences": count,
                "n_rows": rows_containing[smiles],
                "fraction_of_rows": round(rows_containing[smiles] / len(df), 4),
            }
        )
    unique = pd.DataFrame(records).sort_values(
        ["n_occurrences", "exact_mass"], ascending=[False, True]
    ).reset_index(drop=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    unique.to_csv(args.out, index=False)

    total = sum(occurrences.values())
    print(f"Rows read: {len(df)}; fragments: {total}; unique products: {len(unique)}; invalid fragments skipped: {invalid}")
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
