"""Generate a DESL dataset CSV by running the lignin pathway model on .mol files.

Each .mol file under <mol-dir>/<dp>/ is one input lignin structure (folder name
= degree of polymerization). Every file is run through
utils.desl_lignin_pathways.apply_acidolysis and one row is written with the
same columns as data_analysis/df_data_desl.csv.

Sample command:
python data_analysis/generate_desl_dataset.py --dp-min 15 --dp-max 25 \
    --cleavage-fraction 0.9 --alpha-o4-cleavage-fraction 1.0 \
    --acyclic-alpha-o4-cleavage-fraction 1.0 --dbdo-cleavage-fraction 1.0 \
    --out df_data_desl_pathways_dp15_25_cleave90_rings100.csv
"""

import argparse
import sys
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from desl_lignin_pathways import (  # noqa: E402
    apply_acidolysis,
    resolve_route_weights,
)
from desl_simulation import (  # noqa: E402
    calculate_dp,
    calculate_pdi,
    find_functional_groups_counts,
)
from rdkit.Chem import Descriptors  # noqa: E402


def build_parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--mol-dir", default="data_analysis/mol")
    p.add_argument("--dp-min", type=int, default=18)
    p.add_argument("--dp-max", type=int, default=25)
    p.add_argument("--cleavage-fraction", type=float, default=0.8)
    p.add_argument("--p-hibbert", type=float, default=None,
                   help="Hibbert route weight (default 0.57; normalized).")
    p.add_argument("--p-c2-aldehyde", type=float, default=None,
                   help="C2-aldehyde route weight (default 0.40; normalized).")
    p.add_argument("--p-vanillin", type=float, default=None,
                   help="Vanillin-type route weight (default 0.02; normalized).")
    p.add_argument("--vanillin-branch-fraction", type=float, default=None,
                   help="Legacy: vanillin probability, Hibbert otherwise. "
                        "Cannot be combined with --p-* options.")
    p.add_argument("--alpha-o4-cleavage-fraction", type=float, default=0.5)
    p.add_argument("--acyclic-alpha-o4-cleavage-fraction", type=float, default=0.5)
    p.add_argument("--dbdo-cleavage-fraction", type=float, default=0.5)
    p.add_argument("--random-seed", type=int, default=10)
    p.add_argument(
        "--keep-duplicates",
        action="store_true",
        help="Simulate repeated input structures too (default: unique inputs only).",
    )
    p.add_argument("--out", default="df_data_desl_pathways_dp15_25_cleave90_rings100.csv")
    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    try:
        route_weights = resolve_route_weights(
            hibbert=args.p_hibbert,
            c2_aldehyde=args.p_c2_aldehyde,
            vanillin=args.p_vanillin,
            vanillin_branch_fraction=args.vanillin_branch_fraction,
        )
    except ValueError as error:
        parser.error(str(error))
    RDLogger.DisableLog("rdApp.*")
    rows, skipped, seen = [], [], set()
    duplicate_inputs = 0

    for dp in range(args.dp_min, args.dp_max + 1):
        folder = Path(args.mol_dir) / str(dp)
        if not folder.is_dir():
            print(f"Skipping missing folder: {folder}")
            continue
        files = sorted(folder.glob("*.mol"))
        for i, path in enumerate(files, 1):
            try:
                mol = Chem.MolFromMolFile(str(path))
                if mol is None:
                    raise ValueError("unreadable mol file")
                smiles_in = Chem.MolToSmiles(mol)
                key = Chem.MolToSmiles(mol, isomericSmiles=False)
                if not args.keep_duplicates:
                    if key in seen:
                        duplicate_inputs += 1
                        continue
                    seen.add(key)
                products = apply_acidolysis(
                    smiles_in,
                    cleavage_fraction=args.cleavage_fraction,
                    route_weights=route_weights,
                    random_seed=args.random_seed,
                    acyclic_alpha_o4_cleavage_fraction=args.acyclic_alpha_o4_cleavage_fraction,
                    alpha_o4_cleavage_fraction=args.alpha_o4_cleavage_fraction,
                    dbdo_cleavage_fraction=args.dbdo_cleavage_fraction,
                )[0][0]
                smi = ".".join(products)
                wt_avg, num_avg, pdi = calculate_pdi(smi)
                row = {
                    "smiles": smi,
                    "dataset": "DESL",
                    "Mw": int(Descriptors.ExactMolWt(Chem.MolFromSmiles(smi))),
                    "Wt_avg_Mw": wt_avg,
                    "Num_avg_Mw": num_avg,
                    "PDI": pdi,
                }
                row.update(find_functional_groups_counts(smi))
                row["dp_desl"] = calculate_dp(smi)
                row["dp"] = dp
                row["source_file"] = path.name
                rows.append(row)
            except Exception as exc:  # keep going; report at the end
                skipped.append((str(path), f"{type(exc).__name__}: {exc}"))
            if i % 100 == 0:
                print(f"dp={dp}: {i}/{len(files)}", flush=True)

    df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out)
    print(
        f"\nWrote {len(df)} rows to {args.out}; skipped {len(skipped)} files; "
        f"{duplicate_inputs} duplicate input structures ignored; "
        f"{df['smiles'].duplicated().sum()} duplicate product mixtures."
    )
    for path, reason in skipped[:10]:
        print(f"  skipped {path}: {reason}")


if __name__ == "__main__":
    main()
