"""Rule-based lignin linkage-cleavage and product-scenario model using RDKit.

Transformations run in this order:

1. Optional hydrolysis of acyclic alpha-O-4 aryl ethers.
2. Optional hydrolysis of eight-membered DBDO-like cyclic ethers.
3. Optional hydrolysis of five-membered phenylcoumaran cyclic ethers.
4. Beta-O-4 cleavage, with each cleaved site randomly assigned to the
   Hibbert-ketone route, the C2-aldehyde route (formaldehyde loss), or an
   optional vanillin-type pathway, using normalized route weights.

The alpha-O-4 and DBDO ring-opening rules are assumptions based on literature
and are disabled by default. Cleavage fractions and route weights are scenario
inputs, not kinetic predictions; the default weights (0.57 / 0.40 / 0.02) are
assumptions based on literature product patterns, used as per-site
probabilities. The vanillin-type pathway is included as a low-frequency
pathway, an assumption based on literature reports of vanillin as a minor
product; it includes phenol and acetaldehyde coproducts for atom balance.

The optional molecular-weight precipitation cutoff is a coarse partition
model, not a measured solubility model. Formula-based Mn/Mw/PDI and total-OH
values are not direct GPC or 31P NMR predictions.
"""

import argparse
import random

from rdkit import Chem
from rdkit.Chem import Descriptors, rdChemReactions

try:
    from desl_simulation import find_functional_groups_counts
except ImportError:  # imported as part of the utils package
    from utils.desl_simulation import find_functional_groups_counts


DEFAULT_BETA_O4_SMILES = "COc1cc(C(O)C(CO)Oc2ccc(C(O)C(CO)Oc3ccc(C(O)C(CO)Oc4ccc(C(O)C(CO)Oc5ccc(C(O)C(CO)Oc6ccc(C7OCC8C(c9cc(OC)c%10c(c9)-c9cc(C%11OCC%12C(c%13cc(OC)c%14c(c%13)C(CO)C(c%13ccc(OC(CO)C(O)c%15ccc(O)c(OC)c%15)c(OC)c%13)O%14)OCC%11%12)cc(OC)c9OC(CO)C(c9cc(OC)c%11c(c9)C(CO)C(c9ccc(O)c(OC)c9)O%11)O%10)OCC78)cc6OC)cc5OC)cc4OC)cc3OC)cc2OC)ccc1O"
#DEFAULT_BETA_O4_SMILES = "O(Cc1ccccc1)c2ccccc2"

HIBBERT_KETONE_SMARTS = (
    "[c:1]-[CH:2]([OH:3])-[CH:4]([O:5]-[c:6])-[CH2:7][OH:8]"
    ">>"
    "[c:1]-[CH2:2]-[C:4](=[O:3])-[CH2:7][OH:8].[c:6]-[OH:5]"
)
VANILLIN_TYPE_SMARTS = (
    "[c:1]-[CH:2]([OH:3])-[CH:4]([O:5]-[c:6])-[CH2:7][OH:8]"
    ">>"
    "[c:1]-[CH:2]=[O:3].[c:6]-[OH:5].[CH3:7]-[CH:4]=[O:8]"
)
C2_ALDEHYDE_SMARTS = (
    "[c:1]-[CH:2]([OH:3])-[CH:4]([O:5]-[c:6])-[CH2:7][OH:8]"
    ">>"
    "[c:1]-[CH2:2]-[CH:4]=[O:3].[c:6]-[OH:5].[CH2:7]=[O:8]"
)
# Beta-O-4 units without a gamma-CH2OH (e.g. 2-phenoxy-1-phenylethanol) cannot
# follow the three routes above; they give an aryl-acetaldehyde plus phenol
# (the C2-aldehyde route with no formaldehyde to lose).
NO_GAMMA_ALDEHYDE_SMARTS = (
    "[c:1]-[CH:2]([OH:3])-[CH2:4]-[O:5]-[c:6]"
    ">>"
    "[c:1]-[CH2:2]-[CH:4]=[O:3].[c:6]-[OH:5]"
)
# Default beta-O-4 route weights (Hibbert, C2-aldehyde, vanillin-type): assumptions
# based on literature, normalized before use.
DEFAULT_ROUTE_WEIGHTS = (0.57, 0.40, 0.02)
ACYCLIC_ALPHA_O4_QUERY_SMARTS = "[c:1]-[C:2]-[O:3]-[c:4]"
PHENYLCOUMARAN_ALPHA_O4_QUERY_SMARTS = "[c:1]-[CH:2]([O:3]-[c:4])-[C:5]"


def normalize_route_weights(hibbert, c2_aldehyde, vanillin):
    """Return (hibbert, c2_aldehyde, vanillin) scaled to sum to 1.

    Weights must be non-negative and not all zero.
    """
    weights = (hibbert, c2_aldehyde, vanillin)
    if any(weight < 0 for weight in weights):
        raise ValueError("Route weights must be non-negative")
    total = sum(weights)
    if total <= 0:
        raise ValueError("At least one route weight must be greater than zero")
    return tuple(weight / total for weight in weights)


def resolve_route_weights(hibbert=None, c2_aldehyde=None, vanillin=None, vanillin_branch_fraction=None):
    """Resolve CLI-style route settings into normalized weights.

    The legacy vanillin_branch_fraction (Hibbert otherwise, no C2-aldehyde) cannot
    be combined with the explicit route weights. Unset weights take the defaults.
    """
    explicit = [value is not None for value in (hibbert, c2_aldehyde, vanillin)]
    if vanillin_branch_fraction is not None:
        if any(explicit):
            raise ValueError(
                "Use either vanillin_branch_fraction or the explicit route weights, not both"
            )
        if not 0.0 <= vanillin_branch_fraction <= 1.0:
            raise ValueError("vanillin_branch_fraction must be between 0.0 and 1.0")
        return (1.0 - vanillin_branch_fraction, 0.0, vanillin_branch_fraction)
    defaults = DEFAULT_ROUTE_WEIGHTS
    return normalize_route_weights(
        defaults[0] if hibbert is None else hibbert,
        defaults[1] if c2_aldehyde is None else c2_aldehyde,
        defaults[2] if vanillin is None else vanillin,
    )


def build_parser():
    """Build the CLI parser for input structure and assumed model settings."""
    parser = argparse.ArgumentParser(
        description=(
            "Apply optional alpha-O-4 hydrolysis (an assumption based on literature) "
            "followed by beta-O-4 cleavage with selectable Hibbert, C2-aldehyde or "
            "vanillin-type products."
        )
    )
    parser.add_argument(
        "--smiles",
        default=DEFAULT_BETA_O4_SMILES,
        help="Input lignin or model-compound structure as SMILES.",
    )
    parser.add_argument(
        "--cleavage-fraction",
        type=float,
        default=0.8,
        help=(
            "Assumed fraction of beta-O-4 matches to cleave after optional ring "
            "opening (default: 0.8; not a kinetic prediction)."
        ),
    )
    parser.add_argument(
        "--p-hibbert",
        type=float,
        default=None,
        help=(
            "Weight of the Hibbert-ketone route per cleaved beta-O-4 site "
            "(default: 0.57). Route weights are normalized to sum to 1; "
            "assumptions based on literature."
        ),
    )
    parser.add_argument(
        "--p-c2-aldehyde",
        type=float,
        default=None,
        help=(
            "Weight of the C2-aldehyde route with formaldehyde loss "
            "(default: 0.40; assumption based on literature)."
        ),
    )
    parser.add_argument(
        "--p-vanillin",
        type=float,
        default=None,
        help=(
            "Weight of the vanillin-type pathway (default: 0.02; low-frequency "
            "pathway, an assumption based on literature; also gives phenol and "
            "acetaldehyde)."
        ),
    )
    parser.add_argument(
        "--vanillin-branch-fraction",
        type=float,
        default=None,
        help=(
            "Legacy option: probability of the vanillin-type pathway, Hibbert "
            "otherwise, no C2-aldehyde. Cannot be combined with --p-* options."
        ),
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=10,
        help="Seed for reproducible pathway draws (default: 10).",
    )
    parser.add_argument(
        "--alpha-o4-cleavage-fraction",
        type=float,
        default=0.0,
        help=(
            "Fraction of phenylcoumaran alpha-O-4 rings to hydrolyze (assumption based on literature) "
            "(default: 0; targets five-membered rings; consumes water)."
        ),
    )
    parser.add_argument(
        "--acyclic-alpha-o4-cleavage-fraction",
        type=float,
        default=0.0,
        help=(
            "Fraction of non-ring alpha-O-4 aryl ethers to hydrolyze (assumption based on literature) "
            "(default: 0; consumes water)."
        ),
    )
    parser.add_argument(
        "--dbdo-cleavage-fraction",
        type=float,
        default=0.0,
        help=(
            "Fraction of eight-membered DBDO-like alpha-O-4 rings to hydrolyze "
            "(assumption based on literature; default: 0; ring-composition filter; "
            "consumes water)."
        ),
    )
    parser.add_argument(
        "--reaction-time-hours",
        type=float,
        default=2.0,
        help="Reaction time metadata to report; it does not determine conversion.",
    )
    parser.add_argument(
        "--min-precipitate-mw",
        type=float,
        default=None,
        metavar="DA",
        help=(
            "Optional illustrative molecular-weight cutoff in Da; fragments at "
            "or above it are classified as precipitated DESL. Calibrate against "
            "recovery data; this is not a solubility prediction."
        ),
    )
    return parser


def matches_cyclic_alpha_o4(component, query, ring_size):
    """Return alpha-aryl-ether matches with the requested cyclic-ring topology.

    Five-membered matches must have a beta-carbon aromatic connection
    (phenylcoumaran-like); eight-membered matches use the DBDO-like elemental
    ring composition checked below.
    """
    matches = []
    for match in component.GetSubstructMatches(query):
        alpha_carbon_idx, ether_oxygen_idx, _, beta_carbon_idx = match[1:]
        ether_bond = component.GetBondBetweenAtoms(alpha_carbon_idx, ether_oxygen_idx)
        candidate_rings = [
            ring
            for ring in component.GetRingInfo().AtomRings()
            if len(ring) == ring_size
            and alpha_carbon_idx in ring
            and ether_oxygen_idx in ring
        ]
        for ring in candidate_rings:
            ring_atoms = [component.GetAtomWithIdx(index) for index in ring]
            oxygen_count = sum(atom.GetAtomicNum() == 8 for atom in ring_atoms)
            aromatic_carbon_count = sum(
                atom.GetAtomicNum() == 6 and atom.GetIsAromatic()
                for atom in ring_atoms
            )
            aliphatic_carbon_count = sum(
                atom.GetAtomicNum() == 6 and not atom.GetIsAromatic()
                for atom in ring_atoms
            )
            beta_carbon = component.GetAtomWithIdx(beta_carbon_idx)
            beta_is_aryl_connected = any(
                neighbor.GetIsAromatic()
                and neighbor.GetIdx() != alpha_carbon_idx
                for neighbor in beta_carbon.GetNeighbors()
            )

            if ring_size == 5:
                expected_composition = (1, 2, 2)
                topology_matches = beta_is_aryl_connected
            else:
                expected_composition = (2, 4, 2)
                topology_matches = True

            if (
                ether_bond.IsInRing()
                and (oxygen_count, aromatic_carbon_count, aliphatic_carbon_count)
                == expected_composition
                and topology_matches
            ):
                matches.append(match)
                break
    return matches


def apply_cyclic_alpha_o4_hydrolysis(components, cleavage_fraction, ring_size, label):
    """Hydrolyze a fraction of matching cyclic alpha-O-4 ether bonds in place.

    One water molecule is included per event. The C-alpha--O bond is removed,
    water oxygen is attached to C-alpha, and the original ether oxygen remains
    attached to its aromatic carbon. Each event is sanitized and mass-checked.

    Returns the updated component list, initial matching-site count, and number
    of applied events. If the input has no matching ring, the components are
    returned unchanged with zero counts.
    """
    query = Chem.MolFromSmarts(PHENYLCOUMARAN_ALPHA_O4_QUERY_SMARTS)
    if query is None:
        raise RuntimeError(f"Could not compile the {label} alpha-O-4 query")

    site_count = sum(
        len(matches_cyclic_alpha_o4(component, query, ring_size))
        for component in components
    )
    if site_count == 0 and cleavage_fraction > 0.0:
        return components, 0, 0

    cleavage_count = int(site_count * cleavage_fraction + 0.5)
    water = Chem.MolFromSmiles("O")
    water_mass = Descriptors.ExactMolWt(water)
    for _ in range(cleavage_count):
        changed_component = False
        for component_index, component in enumerate(components):
            matches = matches_cyclic_alpha_o4(component, query, ring_size)
            if not matches:
                continue

            _, alpha_carbon_idx, ether_oxygen_idx, _, _ = matches[0]
            editable = Chem.RWMol(component)
            editable.RemoveBond(alpha_carbon_idx, ether_oxygen_idx)
            water_oxygen_idx = editable.AddAtom(Chem.Atom(8))
            editable.AddBond(alpha_carbon_idx, water_oxygen_idx, Chem.BondType.SINGLE)
            product = editable.GetMol()
            Chem.SanitizeMol(product)

            reactant_mass = Descriptors.ExactMolWt(component) + water_mass
            product_mass = Descriptors.ExactMolWt(product)
            if abs(reactant_mass - product_mass) > 1e-6:
                raise ValueError(f"{label} alpha-O-4 hydrolysis failed the mass-balance check")

            components[component_index:component_index + 1] = list(
                Chem.GetMolFrags(product, asMols=True, sanitizeFrags=True)
            )
            changed_component = True
            break

        if not changed_component:
            raise RuntimeError(f"Could not apply the requested {label} ring openings")

    return components, site_count, cleavage_count


def acyclic_alpha_o4_matches(component, query):
    """Return non-ring aryl-C(alpha)-O-aryl matches with alpha CH or CH2."""
    matches = []
    for match in component.GetSubstructMatches(query):
        alpha_carbon_idx, ether_oxygen_idx = match[1], match[2]
        alpha_carbon = component.GetAtomWithIdx(alpha_carbon_idx)
        ether_bond = component.GetBondBetweenAtoms(alpha_carbon_idx, ether_oxygen_idx)
        if (
            not ether_bond.IsInRing()
            and alpha_carbon.GetTotalNumHs() in (1, 2)
        ):
            matches.append(match)
    return matches


def apply_acyclic_alpha_o4_hydrolysis(components, cleavage_fraction):
    """Hydrolyze non-cyclic alpha-O-4 aryl ethers using one water per event.

    The new water oxygen is attached to the alpha carbon, while the original
    ether oxygen becomes phenolic OH. Returns updated components and target and
    applied-site counts. Missing motifs are skipped without error.
    """
    query = Chem.MolFromSmarts(ACYCLIC_ALPHA_O4_QUERY_SMARTS)
    if query is None:
        raise RuntimeError("Could not compile the acyclic alpha-O-4 query")

    site_count = sum(
        len(acyclic_alpha_o4_matches(component, query)) for component in components
    )
    cleavage_count = int(site_count * cleavage_fraction + 0.5)
    water = Chem.MolFromSmiles("O")
    water_mass = Descriptors.ExactMolWt(water)

    for _ in range(cleavage_count):
        changed_component = False
        for component_index, component in enumerate(components):
            matches = acyclic_alpha_o4_matches(component, query)
            if not matches:
                continue

            _, alpha_carbon_idx, ether_oxygen_idx, _ = matches[0]
            editable = Chem.RWMol(component)
            editable.RemoveBond(alpha_carbon_idx, ether_oxygen_idx)
            water_oxygen_idx = editable.AddAtom(Chem.Atom(8))
            editable.AddBond(alpha_carbon_idx, water_oxygen_idx, Chem.BondType.SINGLE)
            product = editable.GetMol()
            Chem.SanitizeMol(product)

            reactant_mass = Descriptors.ExactMolWt(component) + water_mass
            if abs(reactant_mass - Descriptors.ExactMolWt(product)) > 1e-6:
                raise ValueError("Acyclic alpha-O-4 hydrolysis failed the mass-balance check")

            components[component_index:component_index + 1] = list(
                Chem.GetMolFrags(product, asMols=True, sanitizeFrags=True)
            )
            changed_component = True
            break

        if not changed_component:
            raise RuntimeError("Could not apply the requested acyclic alpha-O-4 cleavages")

    return components, site_count, cleavage_count


def apply_acidolysis(
    smiles,
    cleavage_fraction=0.8,
    vanillin_branch_fraction=0.0,
    random_seed=10,
    acyclic_alpha_o4_cleavage_fraction=0.0,
    alpha_o4_cleavage_fraction=0.0,
    dbdo_cleavage_fraction=0.0,
    route_weights=None,
):
    """Apply optional ring openings, then beta-O-4 cleavage, to a SMILES input.

    Processing order is acyclic alpha-O-4 hydrolysis, DBDO-like ring opening,
    phenylcoumaran ring opening, then beta-O-4 cleavage. For each cleaved
    beta-O-4 site, a seeded random draw selects the vanillin-type pathway, the
    C2-aldehyde route (formaldehyde loss), or the Hibbert route using the
    normalized route weights. The vanillin route includes phenol and
    acetaldehyde coproducts (an assumption based on literature). Site selection follows RDKit
    match order, not measured kinetics.

    Args:
        smiles: Input lignin or model-compound SMILES.
        cleavage_fraction: Assumed fraction of post-ring-opening beta-O-4 sites.
        vanillin_branch_fraction: Legacy setting used only when route_weights is
            None: probability of the vanillin-type pathway, Hibbert otherwise.
        route_weights: Optional (hibbert, c2_aldehyde, vanillin) weights, normalized
            internally. Assumptions based on literature.
        random_seed: Seed used to make pathway draws reproducible.
        acyclic_alpha_o4_cleavage_fraction: Assumed fraction of non-ring alpha-O-4
            aryl ethers to hydrolyze. This pathway is an assumption based on literature
            and consumes water.
        alpha_o4_cleavage_fraction: Assumed fraction of matching phenylcoumaran
            rings to hydrolyze. This pathway is an assumption based on literature and
            consumes water.
        dbdo_cleavage_fraction: Assumed fraction of matching DBDO-like rings to
            hydrolyze. This pathway is an assumption based on literature and consumes
            water.

    Returns:
        A tuple containing one tuple of product SMILES and the matched/applied
        counts for beta-O-4 sites, vanillin and Hibbert products, acyclic alpha-O-4
        sites, phenylcoumaran rings, DBDO-like rings, and (appended last) the
        number of C2-aldehyde events, then the matched and applied counts of
        beta-O-4 units without a gamma-CH2OH (aryl-acetaldehyde + phenol).
    """
    substrate = Chem.MolFromSmiles(smiles)
    if substrate is None:
        raise ValueError("Input is not valid SMILES")
    if not 0.0 <= cleavage_fraction <= 1.0:
        raise ValueError("cleavage_fraction must be between 0.0 and 1.0")
    if route_weights is None:
        route_weights = resolve_route_weights(
            vanillin_branch_fraction=vanillin_branch_fraction
        )
    else:
        route_weights = normalize_route_weights(*route_weights)
    hibbert_weight, c2_weight, vanillin_weight = route_weights
    if not 0.0 <= acyclic_alpha_o4_cleavage_fraction <= 1.0:
        raise ValueError(
            "acyclic_alpha_o4_cleavage_fraction must be between 0.0 and 1.0"
        )
    if not 0.0 <= alpha_o4_cleavage_fraction <= 1.0:
        raise ValueError("alpha_o4_cleavage_fraction must be between 0.0 and 1.0")
    if not 0.0 <= dbdo_cleavage_fraction <= 1.0:
        raise ValueError("dbdo_cleavage_fraction must be between 0.0 and 1.0")

    components = list(Chem.GetMolFrags(substrate, asMols=True, sanitizeFrags=True))
    components, acyclic_alpha_site_count, acyclic_alpha_cleavage_count = (
        apply_acyclic_alpha_o4_hydrolysis(
            components, acyclic_alpha_o4_cleavage_fraction
        )
    )
    components, dbdo_site_count, dbdo_cleavage_count = (
        apply_cyclic_alpha_o4_hydrolysis(
            components, dbdo_cleavage_fraction, 8, "DBDO-like"
        )
    )
    components, alpha_site_count, alpha_cleavage_count = (
        apply_cyclic_alpha_o4_hydrolysis(
            components, alpha_o4_cleavage_fraction, 5, "phenylcoumaran"
        )
    )

    target_query = Chem.MolFromSmarts(HIBBERT_KETONE_SMARTS.split(">>")[0])
    if target_query is None:
        raise RuntimeError("Could not compile the beta-O-4 reactant query")
    site_count = sum(
        len(component.GetSubstructMatches(target_query)) for component in components
    )

    cleavage_count = int(site_count * cleavage_fraction + 0.5)
    random_generator = random.Random(random_seed)
    vanillin_count = 0
    hibbert_count = 0
    c2_count = 0
    for _ in range(cleavage_count):
        draw = random_generator.random()
        if draw < vanillin_weight:
            reaction_smarts = VANILLIN_TYPE_SMARTS
            pathway_name = "vanillin-type pathway"
            vanillin_count += 1
        elif draw < vanillin_weight + c2_weight:
            reaction_smarts = C2_ALDEHYDE_SMARTS
            pathway_name = "C2-aldehyde pathway"
            c2_count += 1
        else:
            reaction_smarts = HIBBERT_KETONE_SMARTS
            pathway_name = "Hibbert pathway"
            hibbert_count += 1

        reaction = rdChemReactions.ReactionFromSmarts(reaction_smarts)
        if reaction is None:
            raise RuntimeError(f"Could not compile the {pathway_name} template")
        reaction.Initialize()

        changed_component = False
        for component_index, component in enumerate(components):
            outcomes = reaction.RunReactants((component,), maxProducts=1)
            if not outcomes:
                continue

            reaction_products = list(outcomes[0])
            for product in reaction_products:
                Chem.SanitizeMol(product)
            product_mass = sum(
                Descriptors.ExactMolWt(product) for product in reaction_products
            )
            if abs(Descriptors.ExactMolWt(component) - product_mass) > 1e-6:
                raise ValueError(
                    f"{pathway_name} failed the molecular mass-balance check"
                )

            components[component_index:component_index + 1] = reaction_products
            changed_component = True
            break

        if not changed_component:
            raise RuntimeError(f"Could not apply the requested {pathway_name} event")

    # Beta-O-4 units without a gamma-CH2OH: single route, same cleavage fraction,
    # no random draw (so the draws above are unaffected).
    no_gamma_query = Chem.MolFromSmarts(NO_GAMMA_ALDEHYDE_SMARTS.split(">>")[0])
    no_gamma_reaction = rdChemReactions.ReactionFromSmarts(NO_GAMMA_ALDEHYDE_SMARTS)
    if no_gamma_query is None or no_gamma_reaction is None:
        raise RuntimeError("Could not compile the no-gamma aldehyde template")
    no_gamma_reaction.Initialize()
    no_gamma_site_count = sum(
        len(component.GetSubstructMatches(no_gamma_query)) for component in components
    )
    no_gamma_cleavage_count = int(no_gamma_site_count * cleavage_fraction + 0.5)
    for _ in range(no_gamma_cleavage_count):
        changed_component = False
        for component_index, component in enumerate(components):
            outcomes = no_gamma_reaction.RunReactants((component,), maxProducts=1)
            if not outcomes:
                continue
            reaction_products = list(outcomes[0])
            for product in reaction_products:
                Chem.SanitizeMol(product)
            product_mass = sum(
                Descriptors.ExactMolWt(product) for product in reaction_products
            )
            if abs(Descriptors.ExactMolWt(component) - product_mass) > 1e-6:
                raise ValueError(
                    "no-gamma aldehyde pathway failed the molecular mass-balance check"
                )
            components[component_index:component_index + 1] = reaction_products
            changed_component = True
            break
        if not changed_component:
            raise RuntimeError("Could not apply the requested no-gamma aldehyde event")

    product_smiles = []
    for fragment in components:
        for atom in fragment.GetAtoms():
            atom.SetAtomMapNum(0)
        product_smiles.append(Chem.MolToSmiles(fragment, canonical=True))

    return (
        [tuple(sorted(product_smiles))],
        site_count,
        cleavage_count,
        vanillin_count,
        hibbert_count,
        acyclic_alpha_site_count,
        acyclic_alpha_cleavage_count,
        alpha_site_count,
        alpha_cleavage_count,
        dbdo_site_count,
        dbdo_cleavage_count,
        c2_count,
        no_gamma_site_count,
        no_gamma_cleavage_count,
    )


def summarize_fragments(smiles_list):
    """Summarize a list of product SMILES using molecular formula descriptors.

    Mn, Mw, and PDI are calculated from exact molecular weights. Total OH
    counts all explicit or implicit O-H atoms and is reported as mmol/g; this
    is only a graph-based proxy, not functional-group-resolved 31P NMR data.
    """
    molecular_weights = []
    hydroxyl_count = 0
    for smiles in smiles_list:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"Invalid product SMILES: {smiles}")
        molecular_weights.append(Descriptors.ExactMolWt(molecule))
        hydroxyl_count += sum(
            atom.GetAtomicNum() == 8 and atom.GetTotalNumHs(includeNeighbors=True) > 0
            for atom in molecule.GetAtoms()
        )

    if not molecular_weights:
        return {
            "count": 0,
            "total_mass": 0.0,
            "Mn": None,
            "Mw": None,
            "PDI": None,
            "hydroxyl_count": 0,
            "hydroxyl_mmol_per_g": None,
        }

    total_mass = sum(molecular_weights)
    number_average_mw = total_mass / len(molecular_weights)
    weight_average_mw = sum(weight**2 for weight in molecular_weights) / total_mass
    return {
        "count": len(molecular_weights),
        "total_mass": total_mass,
        "Mn": number_average_mw,
        "Mw": weight_average_mw,
        "PDI": weight_average_mw / number_average_mw,
        "hydroxyl_count": hydroxyl_count,
        "hydroxyl_mmol_per_g": 1000.0 * hydroxyl_count / total_mass,
    }


RESINOL_SMARTS = "c[C]1OC[C]2[C](c)OC[C]12"
ETHER_OXYGEN_SMARTS = "[#6;!$([CH3])][OX2][#6;!$([CH3])]"


def count_linkages(smiles):
    """Count lignin linkage motifs in a SMILES string (a dotted string counts all fragments).

    Returns a dictionary with:
    - beta-O-4: units with a free C-alpha OH and a gamma-CH2OH (the standard template);
    - resinol (beta-beta / gamma-O-alpha / alpha-O-gamma): bicyclic tetrahydrofuran units;
    - phenylcoumaran (beta-5 / alpha-O-4): five-membered ring ether sites;
    - DBDO (5-5 / beta-O-4 / alpha-O-4): eight-membered ring ether sites;
    - ether bonds (total): ether oxygens joining two carbons, excluding methoxyl groups.
    The motifs overlap with the total (ring ethers and beta-O-4 oxygens are also ether
    bonds), so the entries are not additive. These are substructure counts, not
    measured linkage abundances.
    """
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise ValueError("Input is not valid SMILES")
    ring_query = Chem.MolFromSmarts(PHENYLCOUMARAN_ALPHA_O4_QUERY_SMARTS)
    beta_o4_query = Chem.MolFromSmarts(HIBBERT_KETONE_SMARTS.split(">>")[0])
    return {
        "beta-O-4": len(molecule.GetSubstructMatches(beta_o4_query)),
        "Resinol (beta-beta'/gamma-O-alpha/alpha-O-gamma)": len(
            molecule.GetSubstructMatches(Chem.MolFromSmarts(RESINOL_SMARTS))
        ),
        "Phenylcoumaran (beta-5/alpha-O-4)": len(
            matches_cyclic_alpha_o4(molecule, ring_query, 5)
        ),
        "Dibenzodioxocin (DBDO) (5-5/beta-O-4/alpha-O-4)": len(
            matches_cyclic_alpha_o4(molecule, ring_query, 8)
        ),
        "ether bonds (total)": len(
            molecule.GetSubstructMatches(Chem.MolFromSmarts(ETHER_OXYGEN_SMARTS))
        ),
    }


def main():
    """Parse arguments, run selected transformations, and print product metrics."""
    parser = build_parser()
    args = parser.parse_args()
    if not 0.0 <= args.cleavage_fraction <= 1.0:
        parser.error("--cleavage-fraction must be between 0 and 1")
    try:
        route_weights = resolve_route_weights(
            hibbert=args.p_hibbert,
            c2_aldehyde=args.p_c2_aldehyde,
            vanillin=args.p_vanillin,
            vanillin_branch_fraction=args.vanillin_branch_fraction,
        )
    except ValueError as error:
        parser.error(str(error))
    if not 0.0 <= args.acyclic_alpha_o4_cleavage_fraction <= 1.0:
        parser.error("--acyclic-alpha-o4-cleavage-fraction must be between 0 and 1")
    if not 0.0 <= args.alpha_o4_cleavage_fraction <= 1.0:
        parser.error("--alpha-o4-cleavage-fraction must be between 0 and 1")
    if not 0.0 <= args.dbdo_cleavage_fraction <= 1.0:
        parser.error("--dbdo-cleavage-fraction must be between 0 and 1")
    if args.reaction_time_hours <= 0:
        parser.error("--reaction-time-hours must be greater than zero")
    if args.min_precipitate_mw is not None and args.min_precipitate_mw <= 0:
        parser.error("--min-precipitate-mw must be greater than zero")

    substrate = Chem.MolFromSmiles(args.smiles)
    if substrate is None:
        raise ValueError("Input is not valid SMILES")

    (
        products,
        site_count,
        cleavage_count,
        vanillin_count,
        hibbert_count,
        acyclic_alpha_site_count,
        acyclic_alpha_cleavage_count,
        alpha_site_count,
        alpha_cleavage_count,
        dbdo_site_count,
        dbdo_cleavage_count,
        c2_count,
        no_gamma_site_count,
        no_gamma_cleavage_count,
    ) = apply_acidolysis(
        args.smiles,
        cleavage_fraction=args.cleavage_fraction,
        route_weights=route_weights,
        random_seed=args.random_seed,
        acyclic_alpha_o4_cleavage_fraction=args.acyclic_alpha_o4_cleavage_fraction,
        alpha_o4_cleavage_fraction=args.alpha_o4_cleavage_fraction,
        dbdo_cleavage_fraction=args.dbdo_cleavage_fraction,
    )

    print("\n\nDES Lignin : Simulation Model")
    print("=============================\n")
    print("DES solvent: choline chloride:lactic acid, 1:2 molar ratio")
    print("Temperature constraint: below 145 C")
    print(f"Assumed reaction time: {args.reaction_time_hours:g} hours")
    print(f"Pathway random seed: {args.random_seed}")
    if args.cleavage_fraction > 0.0 and site_count == 0 and no_gamma_site_count == 0:
        print("No beta-O-4 motifs remain after optional ring opening; beta-O-4 step skipped.")
    if cleavage_count or no_gamma_cleavage_count:
        print("Modeled event: beta-O-4 cleavage with selected product pathways")
    if (
        not cleavage_count
        and not no_gamma_cleavage_count
        and not acyclic_alpha_cleavage_count
        and not alpha_cleavage_count
        and not dbdo_cleavage_count
    ):
        print("No matching cleavage events were applied.")
    if args.acyclic_alpha_o4_cleavage_fraction > 0.0:
        if acyclic_alpha_site_count:
            actual_acyclic_fraction = (
                acyclic_alpha_cleavage_count / acyclic_alpha_site_count
            )
            print("Additional event: acyclic alpha-O-4 hydrolysis (assumption based on literature)")
            print(
                "Requested acyclic alpha-O-4 cleavage fraction: "
                f"{args.acyclic_alpha_o4_cleavage_fraction:.0%}"
            )
            print(
                f"Acyclic alpha-O-4 sites hydrolyzed: "
                f"{acyclic_alpha_cleavage_count}/{acyclic_alpha_site_count} "
                f"({actual_acyclic_fraction:.1%})"
            )
            print("Acyclic alpha-O-4 hydrolysis is an assumption based on literature; water is consumed.")
        else:
            print("No acyclic alpha-O-4 ether matches; optional cleavage skipped.")
    actual_fraction = cleavage_count / site_count if site_count else 0.0
    print(f"Requested beta-O-4 cleavage fraction: {args.cleavage_fraction:.0%}")
    print(
        f"Integer-rounded sites cleaved: {cleavage_count}/{site_count} "
        f"({actual_fraction:.1%})"
    )
    print(
        "Normalized beta-O-4 route weights (Hibbert / C2-aldehyde / vanillin): "
        f"{route_weights[0]:.2f} / {route_weights[1]:.2f} / {route_weights[2]:.2f} "
        "(assumptions based on literature)"
    )
    print(
        f"Randomly assigned beta-O-4 products: {hibbert_count} Hibbert-type, "
        f"{c2_count} C2-aldehyde, {vanillin_count} vanillin-type"
    )
    if cleavage_count:
        print(
            "Actual route shares of cleaved sites: "
            f"{hibbert_count / cleavage_count:.1%} Hibbert, "
            f"{c2_count / cleavage_count:.1%} C2-aldehyde, "
            f"{vanillin_count / cleavage_count:.1%} vanillin"
        )
    if no_gamma_site_count:
        print(
            "Beta-O-4 units without gamma-CH2OH cleaved to aryl-acetaldehyde + phenol: "
            f"{no_gamma_cleavage_count}/{no_gamma_site_count} (assumed route based on literature)"
        )
    if c2_count:
        print(
            "C2-aldehyde route releases formaldehyde; its ratio is an assumption "
            "based on literature."
        )
    if vanillin_count:
        print(
            "Vanillin-type route is a low-frequency pathway (an assumption based on "
            "literature) and gives phenol plus acetaldehyde coproducts for atom balance."
        )
    print("Optional phenylcoumaran and DBDO-like ring opening are assumptions based on literature.")
    if args.alpha_o4_cleavage_fraction > 0.0:
        if alpha_site_count:
            actual_alpha_fraction = alpha_cleavage_count / alpha_site_count
            print("Additional event: phenylcoumaran alpha-O-4 ring hydrolysis")
            print(f"Requested phenylcoumaran ring-opening fraction: {args.alpha_o4_cleavage_fraction:.0%}")
            print(
                f"Integer-rounded phenylcoumaran rings opened: "
                f"{alpha_cleavage_count}/{alpha_site_count} "
                f"({actual_alpha_fraction:.1%})"
            )
            print("One water molecule is consumed per modeled ring-opening event.")
        else:
            print("No phenylcoumaran rings matched; optional cleavage skipped.")
    if args.dbdo_cleavage_fraction > 0.0:
        if dbdo_site_count:
            actual_dbdo_fraction = dbdo_cleavage_count / dbdo_site_count
            print(f"Requested DBDO-like ring-opening fraction: {args.dbdo_cleavage_fraction:.0%}")
            print(
                f"Integer-rounded DBDO-like rings opened: "
                f"{dbdo_cleavage_count}/{dbdo_site_count} "
                f"({actual_dbdo_fraction:.1%})"
            )
            print("DBDO-like opening is an assumption based on literature; one water is consumed per event.")
        else:
            print("No DBDO-like rings matched; optional cleavage skipped.")
    print("Reaction time is metadata; cleavage fraction is not kinetically predicted.")
    if acyclic_alpha_cleavage_count or alpha_cleavage_count or dbdo_cleavage_count:
        print("Water is consumed by the selected alpha-O-4 hydrolysis steps.")
    else:
        print("Water is not a reactant in the beta-O-4 net-reaction template.")
    print("The solvent's effects and any catalytic role of lactic-acid acidity are not modeled.")
    print("Condensation/repolymerization is not modeled.")
    input_smiles = Chem.MolToSmiles(substrate, canonical=True)
    input_summary = summarize_fragments(
        [Chem.MolToSmiles(m) for m in Chem.GetMolFrags(substrate, asMols=True)]
    )
    print("\nINPUT")
    print(f"Molecular structure (SMILES) : \n{input_smiles}")
    print(f"\nMolecular weight (g/mol)           : {int(Descriptors.ExactMolWt(substrate))}")
    print(f"\nFunctional group                   : {find_functional_groups_counts(input_smiles)}")
    print(f"\nPoly dispersity index (PDI)        : {round(input_summary['PDI'], 1)}")

    print("\nOUTPUT")
    water_count = acyclic_alpha_cleavage_count + alpha_cleavage_count + dbdo_cleavage_count
    water_mass = water_count * Descriptors.ExactMolWt(Chem.MolFromSmiles("O"))
    print(
        f"Water molecules added by hydrolysis steps: {water_count} "
        f"(+{water_mass:.1f} g/mol; input {Descriptors.ExactMolWt(substrate):.1f} "
        f"+ water = {Descriptors.ExactMolWt(substrate) + water_mass:.1f})"
    )
    all_fragments = list(products[0])
    all_summary = summarize_fragments(all_fragments)

    if args.min_precipitate_mw is None:
        print("No precipitation cutoff set; all generated fragments are reported.")
        fragment_groups = [("All simulated fragments", all_fragments)]
    else:
        precipitated = []
        soluble = []
        for fragment_smiles in all_fragments:
            molecule = Chem.MolFromSmiles(fragment_smiles)
            molecular_weight = Descriptors.ExactMolWt(molecule)
            if molecular_weight >= args.min_precipitate_mw:
                precipitated.append(fragment_smiles)
            else:
                soluble.append(fragment_smiles)

        precipitated_summary = summarize_fragments(precipitated)
        soluble_summary = summarize_fragments(soluble)
        recovered_fraction = (
            precipitated_summary["total_mass"] / all_summary["total_mass"]
            if all_summary["total_mass"]
            else 0.0
        )
        print(
            f"Illustrative precipitation cutoff: {args.min_precipitate_mw:g} Da "
            "(assumed value)"
        )
        print(
            f"Modeled recovered mass fraction: {recovered_fraction:.1%} "
            "of all simulated product mass"
        )
        fragment_groups = [
            ("Precipitated DESL fraction", precipitated),
            ("Soluble filtrate fraction", soluble),
        ]

    for group_name, fragment_smiles_list in fragment_groups:
        summary = summarize_fragments(fragment_smiles_list)
        print(f"\n{group_name}")
        print(f"\nNumber of fragments  : {summary['count']}")
        if fragment_smiles_list:
            print("\t" + "\n\t".join(fragment_smiles_list))
        if summary["count"]:
            joined = ".".join(fragment_smiles_list)
            print(
                f"\nTotal Molecular weight (g/mol) of all fragments: "
                f"{int(summary['total_mass'])}"
            )
            print(
                f"\nFunctional group                                 : "
                f"{find_functional_groups_counts(joined)}"
            )
            print(
                f"\nPoly dispersity index (PDI)                      : "
                f"{round(summary['PDI'], 1)}"
            )
            print(
                f"\nMn / Mw (Da), formula-based                      : "
                f"{summary['Mn']:.2f} / {summary['Mw']:.2f}"
            )
            print(
                f"\nTotal OH (mmol/g, not a 31P NMR prediction)      : "
                f"{summary['hydroxyl_mmol_per_g']:.2f}"
            )


if __name__ == "__main__":
    main()