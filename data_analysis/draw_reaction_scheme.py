"""Draw the reactions coded in utils/desl_lignin_pathways.py as a scheme image.

Each reaction is applied to a small model compound with the script's own
apply_acidolysis function, so the drawn products are what the code produces.
The ring-opening reactions add one water. Route weights, fractions and the
model compounds are illustrative (see the module docstring of desl_lignin_pathways.py).

Sample command:
python data_analysis/draw_reaction_scheme.py --out image/reaction_scheme.png
"""

import argparse
import io
import sys
from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, rdChemReactions

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "utils"))
from desl_lignin_pathways import apply_acidolysis  # noqa: E402

RDLogger.DisableLog("rdApp.*")

# (title, model-compound SMILES, apply_acidolysis keyword arguments, uses water)
REACTIONS = [
    (
        "1. Acyclic α-O-4 hydrolysis (assumed, + H2O)",
        "c1ccc(COc2ccccc2)cc1",
        {"cleavage_fraction": 0.0, "acyclic_alpha_o4_cleavage_fraction": 1.0},
        True,
    ),
    (
        "2. DBDO-like ring opening (assumed, + H2O)",
        "CCCc1cc(OC)c2c(c1)-c1cc(CCC)cc(OC)c1OC(CO)C(c1ccc(O)c(OC)c1)O2",
        {"cleavage_fraction": 0.0, "dbdo_cleavage_fraction": 1.0},
        True,
    ),
    (
        "3. Phenylcoumaran ring opening (assumed, + H2O)",
        "COc1cc(C2Oc3c(OC)cc(CCCO)cc3C2CO)ccc1O",
        {"cleavage_fraction": 0.0, "alpha_o4_cleavage_fraction": 1.0},
        True,
    ),
    (
        "4. β-O-4 cleavage, Hibbert route",
        "COc1cc(C(O)C(CO)Oc2ccccc2OC)ccc1O",
        {"cleavage_fraction": 1.0, "route_weights": (1, 0, 0)},
        False,
    ),
    (
        "5. β-O-4 cleavage, C2-aldehyde route (formaldehyde loss)",
        "COc1cc(C(O)C(CO)Oc2ccccc2OC)ccc1O",
        {"cleavage_fraction": 1.0, "route_weights": (0, 1, 0)},
        False,
    ),
    (
        "6. β-O-4 cleavage, vanillin-type route (very low frequency)",
        "COc1cc(C(O)C(CO)Oc2ccccc2OC)ccc1O",
        {"cleavage_fraction": 1.0, "route_weights": (0, 0, 1)},
        False,
    ),
    (
        "7. β-O-4 without γ-CH2OH (applied only to a model compound, larger lignin dataset do not include this scenario)",
        "C1=CC=C(C=C1)C(COC2=CC=CC=C2)O",
        {"cleavage_fraction": 1.0},
        False,
    ),
]


# Names of the model compounds above, in the same order (descriptive where the
# compound has no common name)
MODEL_COMPOUND_NAMES = [
    "Benzyl phenyl ether",
    "DBDO (5-5 / \u03b2-O-4 / \u03b1-O-4) model with propyl side chains, guaiacyl at C\u03b1",
    "Phenylcoumaran (\u03b2-5) model: dihydrodehydrodiconiferyl-alcohol type",
    "Guaiacylglycerol-\u03b2-guaiacyl ether (GBG, \u03b2-O-4 model)",
    "Guaiacylglycerol-\u03b2-guaiacyl ether (GBG, \u03b2-O-4 model)",
    "Guaiacylglycerol-\u03b2-guaiacyl ether (GBG, \u03b2-O-4 model)",
    "2-Phenoxy-1-phenylethanol (\u03b2-O-4 model without \u03b3-CH2OH)",
]


def reaction_label(index):
    """Return 'title' plus the model-compound name for the reaction at index."""
    return f"{REACTIONS[index][0]}\nModel compound: {MODEL_COMPOUND_NAMES[index]}"


def reaction_image(reactant, uses_water, kwargs):
    """Return a PIL image of reactant (+ water) -> products from the script."""
    products = apply_acidolysis(reactant, **kwargs)[0][0]
    if tuple(products) == (Chem.CanonSmiles(reactant),):
        raise RuntimeError(f"No reaction occurred for {reactant}")
    left = reactant + (".O" if uses_water else "")
    right = ".".join(products)
    reaction = rdChemReactions.ReactionFromSmarts(f"{left}>>{right}", useSmiles=True)
    options = Draw.MolDrawOptions()
    options.legendFontSize = 18
    return Draw.ReactionToImage(reaction, subImgSize=(300, 260), drawOptions=options)


def build_figure():
    """Return a matplotlib figure with one row per coded reaction."""
    figure, axes = plt.subplots(len(REACTIONS), 1, figsize=(14, 3.5 * len(REACTIONS)))
    for index, (axis, (title, reactant, kwargs, uses_water)) in enumerate(zip(axes, REACTIONS)):
        image = reaction_image(reactant, uses_water, kwargs)
        buffer = io.BytesIO()
        image.save(buffer, format="PNG")
        buffer.seek(0)
        axis.imshow(mpimg.imread(buffer))
        axis.set_title(reaction_label(index), loc="left", fontsize=12, fontweight="bold")
        axis.axis("off")
    figure.suptitle(
        "Reactions coded in desl_lignin_pathways.py (model compounds)",
        fontsize=12,
    )
    figure.tight_layout(rect=[0, 0, 1, 0.98])
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--out", default="image/reaction_scheme.png")
    args = parser.parse_args()
    figure = build_figure()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.out, dpi=200)
    print(f"Saved reaction scheme to {args.out}")


if __name__ == "__main__":
    main()
