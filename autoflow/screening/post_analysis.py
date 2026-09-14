import json
from pathlib import Path
from typing import Dict, Any, List, Tuple, Set
import numpy as np
from scipy.cluster.hierarchy import fclusterdata
from ase import Atoms
from ase.io import read
from ase.data import covalent_radii
from ase.neighborlist import NeighborList


def build_bond_graph(atoms: Atoms, scale: float = 1.2) -> Set[Tuple[int, int]]:
    """Generates an undirected bond graph based on covalent radii."""
    cutoffs = [covalent_radii[atoms[i].number] * scale for i in range(len(atoms))]
    nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
    nl.update(atoms)

    bonds = set()
    for i in range(len(atoms)):
        neighbors, _ = nl.get_neighbors(i)
        for j in neighbors:
            if j > i:
                bonds.add((i, j))
    return bonds


def is_reactive(atoms: Atoms, n_slab_atoms: int, n_ads_ref: int, initial_bonds: Set[Tuple[int, int]], scale: float = 1.2) -> bool:
    """Checks if the adsorbate structure underwent chemical rearrangement."""
    ads = atoms[n_slab_atoms:]
    if len(ads) != n_ads_ref:
        return True

    current_bonds = build_bond_graph(ads, scale=scale)
    for bond in initial_bonds:
        if bond not in current_bonds:
            return True
    return False


def is_adsorbed(atoms: Atoms, n_slab_atoms: int, threshold: float = 3.0) -> Tuple[bool, float]:
    """Calculates minimum distance between slab and adsorbate atoms."""
    slab_pos = atoms[:n_slab_atoms].positions
    ads_pos = atoms[n_slab_atoms:].positions
    diff = ads_pos[:, None, :] - slab_pos[None, :, :]
    dists = np.linalg.norm(diff, axis=2)
    min_dist = float(np.min(dists))
    return bool(min_dist < threshold), min_dist


def is_duplicate(a1: Atoms, a2: Atoms, n_slab_atoms: int, tol: float = 0.5) -> bool:
    """Evaluates RMSD between adsorbate positions across two structures."""
    ads1 = a1[n_slab_atoms:]
    ads2 = a2[n_slab_atoms:]
    if len(ads1) != len(ads2):
        return False
    diff = ads1.positions - ads2.positions
    rmsd = float(np.sqrt((diff**2).sum() / len(ads1)))
    return rmsd < tol


def geom_disagreement(c: Dict[str, Any], others: List[Dict[str, Any]], n_slab_atoms: int) -> Any:
    """Finds minimum RMSD difference relative to structures from other methods."""
    min_rmsd = float("inf")
    for o in others:
        if c["method"] != o["method"]:
            ads1 = c["atoms"][n_slab_atoms:]
            ads2 = o["atoms"][n_slab_atoms:]
            if len(ads1) != len(ads2):
                continue
            diff = ads1.positions - ads2.positions
            rmsd = float(np.sqrt((diff**2).sum() / len(ads1)))
            min_rmsd = min(min_rmsd, rmsd)

    return min_rmsd if min_rmsd != float("inf") else None


def run_post_analysis(
    screening_dir: Path,
    methods: List[str] = None,
    unadsorbed_thresh: float = 3.0,
    max_energy_win: float = 0.8,
    cluster_cutoff: float = 0.05,
    bond_scale: float = 1.2,
    n_select_dft: int = 5
) -> Dict[str, Any]:
    """Analyzes relaxed configurations, checks reactivity/adsorption, and generates output metadata JSON files."""
    if methods is None:
        methods = ["MACE-FT"]

    # 1. Dependency Resolution & Reference Construction
    surface_json = screening_dir / "surface_atoms.json"
    gas_poscar = screening_dir.parent / "gas" / "POSCAR"

    n_slab_atoms = json.loads(surface_json.read_text())["total_atoms"]
    gas_ads = read(gas_poscar)
    n_ads_ref = len(gas_ads)
    initial_bonds = build_bond_graph(gas_ads, scale=bond_scale)

    # 2. Data Loading
    confs = []
    for d in sorted(screening_dir.glob("conf_*")):
        json_file = d / "ensemble_screen.json"
        if not json_file.exists():
            continue

        results = json.loads(json_file.read_text())
        energy_map = {r["method"]: r["energy"] for r in results}

        for method in methods:
            xyz_file = d / f"relaxed_{method}.xyz"
            if xyz_file.exists() and method in energy_map:
                confs.append({
                    "conf": d.name,
                    "method": method,
                    "energy": float(energy_map[method]),
                    "atoms": read(xyz_file)
                })

    if not confs:
        summary = {"n_loaded": 0, "n_filtered": 0, "methods": {}}
        (screening_dir / "ensemble_screening_summary.json").write_text(json.dumps(summary, indent=2))
        return summary

    # 3. Reactivity/Adsorption Evaluation & Energy Window Filtering
    filtered = []
    for method in methods:
        method_confs = [c for c in confs if c["method"] == method]
        if not method_confs:
            continue

        for c in method_confs:
            c["reactive"] = is_reactive(c["atoms"], n_slab_atoms, n_ads_ref, initial_bonds, scale=bond_scale)
            c["adsorbed"], c["ads_dist"] = is_adsorbed(c["atoms"], n_slab_atoms, threshold=unadsorbed_thresh)

        nonreactive = [c for c in method_confs if not c["reactive"]]
        reactive = [c for c in method_confs if c["reactive"]]

        if nonreactive:
            emin_nonreact = min(c["energy"] for c in nonreactive)
            filtered.extend([c for c in nonreactive if c["energy"] - emin_nonreact <= max_energy_win])

        if reactive:
            emin_react = min(c["energy"] for c in reactive)
            filtered.extend([c for c in reactive if c["energy"] - emin_react <= max_energy_win])

    if not filtered:
        summary = {"n_loaded": len(confs), "n_filtered": 0, "methods": {}}
        (screening_dir / "ensemble_screening_summary.json").write_text(json.dumps(summary, indent=2))
        return summary

    # 4. Hierarchical Energy-based Clustering
    clusters = {}
    cluster_id_counter = 1

    for method in methods:
        method_filtered = [c for c in filtered if c["method"] == method]
        if not method_filtered:
            continue

        features = np.array([[c["energy"]] for c in method_filtered])
        names = [f"{c['conf']}|{c['method']}" for c in method_filtered]

        if len(features) == 1:
            labels = np.array([1])
        else:
            labels = fclusterdata(features, t=cluster_cutoff, criterion="distance")

        for name, label in zip(names, labels):
            global_id = int(cluster_id_counter + int(label) - 1)
            clusters.setdefault(global_id, []).append(name)

        cluster_id_counter += int(labels.max())

    # Organize cluster structures and select representatives per method
    method_clusters = {}
    method_reps = {}

    for cluster_id, members in clusters.items():
        _, method = members[0].split("|")
        method_clusters.setdefault(method, {})
        method_clusters[method][str(cluster_id)] = members

    for method, clusts in method_clusters.items():
        reps = []
        for cluster_id, members in clusts.items():
            best = min(
                members,
                key=lambda m: next(
                    c["energy"] for c in filtered
                    if f"{c['conf']}|{c['method']}" == m
                )
            )
            conf_name, method_name = best.split("|")
            entry = next(c for c in filtered if c["conf"] == conf_name and c["method"] == method_name)

            reps.append({
                "cluster": str(cluster_id),
                "conf": conf_name,
                "energy": entry["energy"],
                "reactive": entry["reactive"],
                "adsorbed": entry["adsorbed"]
            })
        method_reps[method] = reps

    # 5. DFT Candidate Selection Pool Assembly
    candidates = []
    for method, reps in method_reps.items():
        for r in reps:
            entry = next(c for c in filtered if c["conf"] == r["conf"] and c["method"] == method)
            candidates.append({
                "conf": r["conf"],
                "method": method,
                "energy": r["energy"],
                "atoms": entry["atoms"],
                "cluster": r["cluster"],
                "reactive": entry["reactive"],
                "adsorbed": entry["adsorbed"],
                "ads_dist": entry["ads_dist"]
            })

    adsorbed_nonreactive = [c for c in candidates if (not c["reactive"]) and c["adsorbed"]]
    unadsorbed_nonreactive = [c for c in candidates if (not c["reactive"]) and not c["adsorbed"]]

    if len(adsorbed_nonreactive) >= 2:
        working_pool = adsorbed_nonreactive + unadsorbed_nonreactive
    elif len(unadsorbed_nonreactive) >= 2:
        working_pool = unadsorbed_nonreactive
    else:
        working_pool = candidates

    # RMSD Deduplication
    unique_candidates = []
    for c in working_pool:
        if not any(is_duplicate(c["atoms"], u["atoms"], n_slab_atoms) for u in unique_candidates):
            unique_candidates.append(c)

    for c in unique_candidates:
        c["disagreement"] = geom_disagreement(c, unique_candidates, n_slab_atoms)

    # Diversity and Enthalpy Priority Selection for DFT
    by_energy = sorted(unique_candidates, key=lambda x: x["energy"])
    by_disagreement = sorted(
        unique_candidates,
        key=lambda x: -x["disagreement"] if x["disagreement"] is not None else float("-inf")
    )

    selected = []

    # Pick lowest energy adsorbed configurations
    for c in by_energy:
        if c["ads_dist"] < unadsorbed_thresh and c not in selected:
            selected.append(c)
        if len(selected) >= 2:
            break

    # Fallback lowest energy
    for c in by_energy:
        if c not in selected:
            selected.append(c)
        if len(selected) >= 2:
            break

    # High geometric disagreement solutions
    has_any_pair = any(c["disagreement"] is not None for c in unique_candidates)
    for c in by_disagreement:
        if c not in selected:
            selected.append(c)
        if len(selected) >= 4:
            break

    # Fallback routine if single method ensemble
    if not has_any_pair:
        methods_present = set(c["method"] for c in selected)
        for method in methods:
            if method not in methods_present:
                pool = [
                    c for c in candidates
                    if c["method"] == method and (not c["reactive"]) and c["ads_dist"] >= unadsorbed_thresh
                ]
                for c in pool:
                    if all(not is_duplicate(c["atoms"], s["atoms"], n_slab_atoms) for s in selected):
                        selected.append(c)
                        break

        mp_pool = [
            c for c in candidates
            if c["method"] == "MACE-MP" and (not c["reactive"]) and c["ads_dist"] < unadsorbed_thresh
        ]
        added = 0
        for c in mp_pool:
            if all(not is_duplicate(c["atoms"], s["atoms"], n_slab_atoms) for s in selected):
                selected.append(c)
                added += 1
            if added >= 2:
                break

    # Fill remaining capacity
    for c in unique_candidates:
        if c not in selected:
            selected.append(c)
        if len(selected) >= n_select_dft:
            break

    # Write dft_selection.json
    dft_output = {
        "n_candidates": len(unique_candidates),
        "n_selected": len(selected),
        "selection": []
    }
    for c in selected:
        item = {
            "conf": c["conf"],
            "method": c["method"],
            "energy": c["energy"],
            "adsorption_distance": c["ads_dist"],
            "reactive": c["reactive"]
        }
        if c.get("disagreement") is not None:
            item["geom_disagreement"] = c["disagreement"]
        dft_output["selection"].append(item)

    (screening_dir / "dft_selection.json").write_text(json.dumps(dft_output, indent=2))

    # 6. Best Candidates Summary Generation
    best_per_method = []
    summary = {
        "n_loaded": len(confs),
        "n_filtered": len(filtered),
        "methods": {}
    }

    for method in methods:
        method_entries = [c for c in filtered if c["method"] == method and c["adsorbed"]]
        if not method_entries:
            continue

        nonreactive = [c for c in method_entries if not c["reactive"]]
        pool = nonreactive if nonreactive else method_entries
        best = min(pool, key=lambda x: x["energy"])

        best_entry = {
            "conf": best["conf"],
            "method": best["method"],
            "energy": best["energy"],
            "reactive": best["reactive"],
            "adsorbed": best["adsorbed"]
        }
        best_per_method.append(best_entry)

        summary["methods"][method] = {
            "n_filtered": len(method_entries),
            "n_clusters": len(method_clusters.get(method, {})),
            "clusters": method_clusters.get(method, {}),
            "representatives": method_reps.get(method, []),
            "best": best_entry
        }

    (screening_dir / "ensemble_screening_summary.json").write_text(json.dumps(summary, indent=2))
    return summary
