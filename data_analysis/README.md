# Data analysis and datasets

This folder holds the dataset generation script, the generated datasets, the analysis notebook and its figures. Run the commands below from the repository root.

- `generate_desl_dataset.py`: runs the model on the LGS structures and writes a dataset CSV.
- `extract_unique_products.py`: lists the unique products in a dataset CSV.
- `draw_reaction_scheme.py`: draws the reaction scheme.
- [datasets/](datasets/): the generated datasets.
- [charts/](charts/): figures and interactive maps made by the analysis notebook.

## Generating a dataset

`data_analysis/generate_desl_dataset.py` runs the model on lignin structures from the LGS dataset and writes one row per structure. **The LGS structures are not included in this repository.** Obtain the LGS dataset from its source ([add the LGS dataset reference or link here]) and place its `.mol` files in `mol/<dp>/` (that is, `data_analysis/mol/<dp>/`), one folder per degree of polymerization (DP). The script then reads the folders from `--dp-min` to `--dp-max`:

```bash
python data_analysis/generate_desl_dataset.py --dp-min 18 --dp-max 25 \
  --cleavage-fraction 1.0 --alpha-o4-cleavage-fraction 1.0 \
  --acyclic-alpha-o4-cleavage-fraction 1.0 --dbdo-cleavage-fraction 1.0 \
  --out data_analysis/datasets/df_data_desl_pathways_100.csv
```

Its defaults are cleavage 0.8, 0.5 for the three optional hydrolyses and the default route weights (0.57 / 0.40 / 0.02, see the [main README](../README.md)). Repeated input structures are skipped unless `--keep-duplicates` is given. All files use the same random seed unless it is changed.

## Generated datasets

The datasets generated with the model are in [datasets/](datasets/). They cover the LGS structures with DP 15-25 (4,411 structures) and use 0.9 β-O-4 cleavage, 1.0 acyclic α-O-4, DBDO and phenylcoumaran hydrolyses, and the default route weights (0.57 / 0.40 / 0.02) and random seed. See [datasets/README.md](datasets/README.md) for full details.

| File | Rows | Description |
|---|---|---|
| [df_data_desl_pathways_dp15_25_cleave90_rings100.csv](datasets/df_data_desl_pathways_dp15_25_cleave90_rings100.csv) | 4,411 | One row per LGS structure with the simulated product mixture (`smiles`), formula-based masses and PDI, and functional-group counts |
| [df_data_desl_pathways_dp15_25_cleave90_rings100_cutoff350.csv](datasets/df_data_desl_pathways_dp15_25_cleave90_rings100_cutoff350.csv) | 4,411 | The same mixtures keeping only fragments of at least 350 g/mol, the approximate precipitated DESL fraction (made with `apply_mass_cutoff.py`) |
| [unique_products_dp15_25_cleave90_rings100.csv](datasets/unique_products_dp15_25_cleave90_rings100.csv) | 482 | The unique products found in the first file (made with `extract_unique_products.py`) |

The analysis notebook [desl_analysis.ipynb](desl_analysis.ipynb) reads these files and saves its figures in [charts/](charts/).
