# Simulated DESL dataset from LGS lignin structures

This dataset contains simulated product mixtures of deep eutectic solvent lignin (DESL), predicted from lignin structures in the LGS dataset. It has three data files:

- `df_data_desl_pathways_dp15_25_cleave90_rings100.csv`: one row per input lignin structure, with the predicted fragment mixture after acidolysis, its masses and dispersity, and functional-group counts.
- `df_data_desl_pathways_dp15_25_cleave90_rings100_cutoff350.csv`: the same 4,411 mixtures, keeping only fragments of at least 350 g/mol. This approximates the precipitated DESL fraction. Masses, dispersity and functional-group counts are calculated from the retained fragments.
- `unique_products_dp15_25_cleave90_rings100.csv`: the 482 unique products found across the 4,411 mixtures in `df_data_desl_pathways_dp15_25_cleave90_rings100.csv`, with formula, exact mass, number of occurrences, and share of mixtures containing each one.

All values are simulated, not measured.

## Citation

> Eswaran, Sudha Cheranma Devi; Zhang, Xiao; Subramaniam, Senthil (2026). *Simulated DESL dataset from LGS lignin structures*. figshare. Dataset. https://doi.org/10.6084/m9.figshare.35113225.v2

## Summary

| Field | Value |
|---|---|
| Title | Simulated DESL dataset from LGS lignin structures |
| Version | 2 |
| DOI | [10.6084/m9.figshare.35113225.v2](https://doi.org/10.6084/m9.figshare.35113225.v2) |
| Repository | figshare |
| Creators | Eswaran, Sudha Cheranma Devi; Zhang, Xiao; Subramaniam, Senthil |
| License | [to be chosen, e.g. CC BY 4.0] |
| Resource type | Dataset, simulated (in silico) |
| Format | CSV, UTF-8, comma-separated, header row |
| Input data | LGS lignin structures, DP 15–25 ([add LGS reference / DOI]) |
| Keywords | lignin, deep eutectic solvent, DESL, acidolysis, β-O-4 cleavage, simulation, SMILES, polydispersity |

## Files

| File | Rows | DP | β-O-4 cleavage | Acyclic α-O-4, DBDO and phenylcoumaran hydrolysis | Description |
|---|---|---|---|---|---|
| `df_data_desl_pathways_dp15_25_cleave90_rings100.csv` | 4,411 | 15–25 | 0.9 | 1.0 | Simulated mixtures, 90% β-O-4 cleavage, all three ether types fully opened |
| `df_data_desl_pathways_dp15_25_cleave90_rings100_cutoff350.csv` | 4,411 | 15–25 | 0.9 | 1.0 | Same mixtures as the `cleave90_rings100` file, keeping only fragments of at least 350 g/mol (approximate precipitated DESL fraction); all columns recalculated from the retained fragments |
| `unique_products_dp15_25_cleave90_rings100.csv` | 482 | 15–25 | 0.9 | 1.0 | Unique products found in `df_data_desl_pathways_dp15_25_cleave90_rings100.csv` |
| `SHA256SUMS` | | | | | SHA-256 checksums of the CSV files |

Rows per DP (both mixture files):

| DP | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 | 25 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Rows | 257 | 296 | 296 | 263 | 311 | 315 | 311 | 315 | 681 | 685 | 681 |

## Simulation conditions

- **Input structures:** Each unique LGS structure was simulated once. Repeated input structures were excluded.
- **Ether hydrolysis:** Before β-O-4 cleavage, acyclic α-O-4 ethers, dibenzodioxocin (DBDO) rings and phenylcoumaran (β-5/α-O-4) rings were hydrolysed at the fractions in the [Files](#files) table.
- **β-O-4 cleavage:** The stated fraction of β-O-4 linkages was cleaved. Each cleavage followed one of three routes, with these relative weights:
  - Hibbert ketone (0.57)
  - C2-aldehyde with formaldehyde release (0.40)
  - vanillin-type with acetaldehyde release (0.02)

  β-O-4 units without a γ-CH₂OH group were cleaved at the same fraction through a single aldehyde-forming route.
- **Route selection:** The number of cleaved linkages is the stated fraction of the available sites, rounded to the nearest integer. Sites are taken in a fixed order, and the route for each cleavage is drawn at random with a fixed seed (10), so the files are reproducible.
- **Mass cutoff (cutoff-350 file only):** Fragments below 350 g/mol were removed from each mixture to approximate the DESL fraction recovered by anti-solvent precipitation. This removes the volatiles (formaldehyde, acetaldehyde), the aromatic monomers and small dimers, about 45% of the product mass on average. The cutoff is a modeling assumption; precipitation is not simulated.
- **Properties:** Masses, dispersity and functional-group counts were calculated from the product SMILES with RDKit (version 2026.03.1).

## Columns: simulated mixture files (`df_data_desl_pathways_*.csv`)

| Column | Type | Unit | Description |
|---|---|---|---|
| (unnamed first column) | integer | | Row index |
| `smiles` | string | | Simulated product mixture: the SMILES of each fragment, joined with `.`. Repeated fragments appear once per copy |
| `dataset` | string | | Always `DESL` |
| `Mw` | integer | g/mol (Da) | Total monoisotopic exact mass of the whole mixture (sum of all fragments), truncated to an integer. This is **not** the weight-average molecular weight |
| `Wt_avg_Mw` | float | g/mol (Da) | Weight-average molecular mass (Mw) of the fragments in the mixture |
| `Num_avg_Mw` | float | g/mol (Da) | Number-average molecular mass (Mn) of the fragments in the mixture |
| `PDI` | float | | Dispersity, `Wt_avg_Mw / Num_avg_Mw` |
| `ether` | integer | count | Matches of SMARTS `[C;H1][OX2][c]` (aryl ether at a CH carbon, e.g. β-O-4) |
| `carbonyl` | integer | count | Matches of SMARTS `[CX3]=[OX1]` |
| `methoxy` | integer | count | Matches of SMARTS `[c][OX2][CH3]` |
| `phenolicOH` | integer | count | Matches of SMARTS `c1ccccc1[OH]` |
| `aliphaticOH_secondary` | integer | count | Matches of SMARTS `cC[OH]` (α-hydroxyl) |
| `aliphaticOH_primary` | integer | count | Matches of SMARTS `cCCC[OH]` (γ-hydroxyl) |
| `cinnamyl_alcohol` | integer | count | Matches of SMARTS `c1ccccc1C=CCO` |
| `dp_desl` | string (list) | | List of the distinct estimated DPs of the fragments, each `int(fragment exact mass / 180)`, e.g. `[0, 1, 2]` |
| `dp` | integer | | DP of the input LGS structure |
| `source_file` | string | | Identifier (file name) of the input LGS structure |

All masses are monoisotopic exact masses, not average molecular weights. Functional-group counts are the number of SMARTS matches across the whole mixture. In the cutoff-350 file, every column (including `Mw`, the total mass) is calculated from the retained fragments only; `dp`, `source_file` and `dataset` are unchanged.

## Columns: `unique_products_dp15_25_cleave90_rings100.csv`

| Column | Type | Unit | Description |
|---|---|---|---|
| `product_smiles` | string | | Canonical SMILES of one distinct product |
| `formula` | string | | Molecular formula |
| `exact_mass` | float | g/mol (Da) | Monoisotopic exact mass |
| `n_occurrences` | integer | count | Total copies of the product across all rows |
| `n_rows` | integer | count | Number of rows (mixtures) that contain the product |
| `fraction_of_rows` | float | | `n_rows / 4,411` |

## Limitations

- Only the reactions listed under [Simulation conditions](#simulation-conditions) are included, at fixed route weights and fractions. Kinetics, condensation, solvent incorporation and phase separation of the products are not modeled.
- The cleavage routes are drawn at random, so each row is one possible outcome for that structure.
- The input structures come from the LGS dataset and share its assumptions.
- The 350 g/mol cutoff is an illustrative choice; the real precipitated fraction depends on the anti-solvent and washing conditions.
- Some input structures give identical product mixtures. Count of repeated `smiles` values: 221 in `df_data_desl_pathways_dp15_25_cleave90_rings100.csv` and 1,634 in the cutoff-350 file, where removing the small fragments makes more mixtures identical.

## File integrity

`SHA256SUMS` lists the SHA-256 checksum of each CSV file.
