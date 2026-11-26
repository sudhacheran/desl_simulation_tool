"""Keep only product fragments at or above a mass cutoff in a DESL dataset CSV.

Each row's `smiles` mixture is filtered to fragments with exact mass >= the
cutoff (an illustrative precipitation cutoff, an assumption). Mw, Wt_avg_Mw,
Num_avg_Mw, PDI, functional-group counts and dp_desl are recalculated from the
retained fragments; dataset, dp and source_file are kept. Rows with no retained
fragment are dropped and reported.

Sample command:
python data_analysis/apply_mass_cutoff.py \
    --in data_analysis/datasets/df_data_desl_pathways_dp15_25_cleave90_rings100.csv \
    --cutoff 350 \
    --out data_analysis/datasets/df_data_desl_pathways_dp15_25_cleave90_rings100_cutoff350.csv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from desl_simulation import (  # noqa: E402
    calculate_dp,
    calculate_pdi,
    find_functional_groups_counts,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--in", dest="infile", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--cutoff", type=float, default=350.0, help="Minimum fragment exact mass (Da).")
    args = parser.parse_args()
    RDLogger.DisableLog("rdApp.*")

    df = pd.read_csv(args.infile, index_col=0)
    mass_cache = {}

    def mass(fragment):
        if fragment not in mass_cache:
            mass_cache[fragment] = Descriptors.ExactMolWt(Chem.MolFromSmiles(fragment))
        return mass_cache[fragment]

    rows, dropped = [], 0
    for _, source in df.iterrows():
        kept = [f for f in source["smiles"].split(".") if mass(f) >= args.cutoff]
        if not kept:
            dropped += 1
            continue
        smi = ".".join(kept)
        wt_avg, num_avg, pdi = calculate_pdi(smi)
        row = {
            "smiles": smi,
            "dataset": source["dataset"],
            "Mw": int(Descriptors.ExactMolWt(Chem.MolFromSmiles(smi))),
            "Wt_avg_Mw": wt_avg,
            "Num_avg_Mw": num_avg,
            "PDI": pdi,
        }
        row.update(find_functional_groups_counts(smi))
        row["dp_desl"] = calculate_dp(smi)
        row["dp"] = source["dp"]
        row["source_file"] = source["source_file"]
        rows.append(row)

    out = pd.DataFrame(rows, columns=df.columns)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out)
    print(f"Wrote {len(out)} rows to {args.out}; dropped {dropped} rows with no fragment >= {args.cutoff:g} Da.")


if __name__ == "__main__":
    main()
