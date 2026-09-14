from pathlib import Path
from typing import List
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from ase.atoms import Atoms
from ase.build import sort
from ase.io import write

SMARTS_PATTERNS = [
    ("alkene", Chem.MolFromSmarts("C=C")),
    ("alkyne", Chem.MolFromSmarts("C#C")),
    ("carbonyl", Chem.MolFromSmarts("[CX3]=[OX1]")),
    ("alcohol", Chem.MolFromSmarts("[OX2H]")),
    ("amine", Chem.MolFromSmarts("[NX3;H2,H1;!$(NC=O)]")),
    ("nitrile", Chem.MolFromSmarts("C#N")),
    ("aromatic_ring", Chem.MolFromSmarts("a1aaaaa1")),
    ("co2", Chem.MolFromSmarts("O=C=O")),
]

def detect_anchor_groups(mol, mol_no_dummy, old_to_new_idx, adsorption_sites):
    """
    Returns anchor groups as lists of atom indices (in mol_no_dummy indexing).

    Automatic selection of point of adsorption is done via filters.
    Three tiers of filters are used:
    1) Open-shell species / intermediates
    2) Closed-shell species with key functional group(s)
    3) Gasteiger charge priority fallback
    """
    anchor_groups = []
    if adsorption_sites:
        for heavy_idx in adsorption_sites:
            anchor_groups.append([old_to_new_idx[heavy_idx]])
        return anchor_groups

    for _, pattern in SMARTS_PATTERNS:
        matches = mol_no_dummy.GetSubstructMatches(pattern)
        for match in matches:
            anchor_groups.append(list(match))

    if anchor_groups:
        return anchor_groups

    try:
        AllChem.ComputeGasteigerCharges(mol_no_dummy)
        charges = [float(a.GetProp("_GasteigerCharge")) for a in mol_no_dummy.GetAtoms()]
        anchor_groups.append([int(np.argmax(np.abs(charges)))])
        return anchor_groups
    except Exception:
        return []

def generate_gas_phase(smiles: str, output_dir: Path) -> List[str]:
    """Generates POSCAR and meta.xyz for a given SMILES string, returning unique element symbols."""
    lg = RDLogger.logger()
    lg.setLevel(RDLogger.CRITICAL)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    mol = Chem.AddHs(mol)
    dummy_atoms = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]

    ads_site_map = {}
    for d in dummy_atoms:
        neighbors = [n for n in d.GetNeighbors() if n.GetAtomicNum() != 0]
        if len(neighbors) != 1:
            raise ValueError("Each dummy atom must bond to exactly one heavy atom")
        ads_site_map.setdefault(neighbors[0].GetIdx(), []).append(d.GetIdx())

    adsorption_sites = list(ads_site_map.keys())

    if len(adsorption_sites) >= 2:
        rw = Chem.RWMol(mol)
        rep_dummies = [ads_site_map[h][0] for h in adsorption_sites]
        for i in range(len(rep_dummies) - 1):
            a1, a2 = rep_dummies[i], rep_dummies[i + 1]
            if not rw.GetBondBetweenAtoms(a1, a2):
                rw.AddBond(a1, a2, Chem.BondType.SINGLE)
        mol = rw.GetMol()

    rw = Chem.RWMol(mol)
    for idx in sorted([a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0], reverse=True):
        rw.RemoveAtom(idx)

    mol_no_dummy = rw.GetMol()
    params = AllChem.ETKDGv3()
    if AllChem.EmbedMolecule(mol_no_dummy, params) != 0:
        raise ValueError("3D embedding of molecule failed")

    conf = mol_no_dummy.GetConformer()
    atoms = Atoms(symbols=[a.GetSymbol() for a in mol_no_dummy.GetAtoms()], positions=conf.GetPositions())

    old_to_new_idx = {}
    j = 0
    for i, a in enumerate(mol.GetAtoms()):
        if a.GetAtomicNum() != 0:
            old_to_new_idx[i] = j
            j += 1

    anchor_groups = detect_anchor_groups(mol, mol_no_dummy, old_to_new_idx, adsorption_sites)
    tags = np.full(len(atoms), -1, dtype=int)
    for group_id, group in enumerate(anchor_groups):
        for idx in group:
            tags[idx] = group_id

    atoms.set_tags(tags)
    sorted_atoms = sort(atoms)
    sorted_atoms.set_pbc(False)
    sorted_atoms.set_cell([20.0, 20.0, 20.0])
    sorted_atoms.center()

    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / "meta.xyz", sorted_atoms)
    write(output_dir / "POSCAR", sorted_atoms, format="vasp")

    # Extract distinct non-dummy chemical elements present in periodic order
    elements = sorted(
        {atom.GetSymbol() for atom in mol_no_dummy.GetAtoms() if atom.GetAtomicNum() > 0},
        key=lambda el: Chem.GetPeriodicTable().GetAtomicNumber(el)
    )
    return elements
