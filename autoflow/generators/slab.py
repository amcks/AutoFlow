import json
from pathlib import Path
from typing import Tuple, Optional, Dict, Any, List
import numpy as np
import spglib
from scipy.spatial import Delaunay
from scipy.cluster.hierarchy import linkage, fcluster
from ase.build import bulk, surface
from ase.constraints import FixAtoms
from ase.neighborlist import NeighborList
from ase.io import read, write
from matplotlib.path import Path as MPLPath


# Helper Functions
def crosses_pbc(indices, frac_coords, threshold=0.8) -> bool:
    """Check if atomic site indices cross periodic boundary conditions."""
    if isinstance(indices, int):
        return False
    coords = frac_coords[np.array(indices), :2]
    for dim in range(2):
        delta = coords[:, dim][:, None] - coords[:, dim][None, :]
        delta -= np.round(delta)
        if np.max(np.abs(delta)) > threshold:
            return True
    return False


def project_onto_plane(points: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Project 3D points onto a plane defined by its normal vector."""
    return points - np.outer(points @ normal, normal)


def is_connected_triangle(tri_atoms: Tuple[int, ...], nl: NeighborList) -> bool:
    """Verify that all atoms in a triangle share mutual neighbor links."""
    for i in tri_atoms:
        neighbors = set(nl.get_neighbors(i)[0])
        if len(neighbors & set(tri_atoms)) < 2:
            return False
    return True


def triangle_max_edge(tri_atoms: Tuple[int, ...], positions: np.ndarray, max_edge: Optional[float] = None) -> bool:
    """Verify the maximum edge length of a triangle candidate."""
    coords = positions[list(tri_atoms)]
    dists = [
        np.linalg.norm(coords[i] - coords[j])
        for i in range(3) for j in range(i + 1, 3)
    ]
    if max_edge is None:
        return True
    return max(dists) <= max_edge


def classify_fcc_hcp(tri_atoms: Tuple[int, ...], positions: np.ndarray, sub_positions: np.ndarray, xy_tol: float = 0.35) -> str:
    """Classify threefold site as HCP or FCC based on subsurface atom proximity."""
    if len(sub_positions) == 0:
        return "fcc"  # Fallback default if single-layer terrace or missing subsurface
    tri_xy = positions[list(tri_atoms), :2]
    centroid_xy = tri_xy.mean(axis=0)
    sub_xy = sub_positions[:, :2]
    distances = np.linalg.norm(sub_xy - centroid_xy, axis=1)
    return "hcp" if np.min(distances) < xy_tol else "fcc"


def cluster_sites(candidates: List[Tuple[int, ...]], positions: np.ndarray, tol: float = 0.5) -> List[Tuple[int, ...]]:
    """Cluster duplicate or near-identical threefold candidate sites."""
    clustered = []
    for cand in candidates:
        coords = positions[list(cand), :2]
        if not clustered:
            clustered.append(cand)
        else:
            dists = [np.linalg.norm(coords.mean(axis=0) - positions[list(c), :2].mean(axis=0)) for c in clustered]
            if all(dd > tol for dd in dists):
                clustered.append(cand)
    return clustered


def build_terrace_adjacency(terrace: List[int], positions: np.ndarray) -> Dict[int, List[int]]:
    """Build planar adjacency graph for terrace atoms."""
    coords = positions[terrace]
    coords_2d = coords[:, :2]

    dists = []
    for i in range(len(coords_2d)):
        for j in range(i + 1, len(coords_2d)):
            d = np.linalg.norm(coords_2d[i] - coords_2d[j])
            dists.append(d)

    if not dists:
        return {}

    dists = np.sort(dists)
    unique = np.unique(np.round(dists, 3))
    if len(unique) < 2:
        return {}

    long = unique[1]
    edge_tol = 1.2 * long

    adjacency = {}
    for idx_i, i in enumerate(terrace):
        adjacency[i] = []
        for idx_j, j in enumerate(terrace):
            if j == i:
                continue
            d = np.linalg.norm(coords_2d[idx_i] - coords_2d[idx_j])
            if d <= edge_tol:
                adjacency[i].append(j)

    return adjacency


def is_rectangular_ring(quad: Tuple[int, ...], positions: np.ndarray, heights: np.ndarray, layer_tol: float, angle_tol: float = 15.0) -> bool:
    """Verify if 4-atom site forms a rectangular geometry on the surface layer."""
    coords = positions[list(quad)]

    zs = [heights[i] for i in quad]
    if max(zs) - min(zs) > layer_tol:
        return False

    center = coords.mean(axis=0)
    vecs = coords - center
    angles = np.arctan2(vecs[:, 1], vecs[:, 0])
    order = np.argsort(angles)
    coords = coords[order]

    edges = [coords[(i + 1) % 4] - coords[i] for i in range(4)]

    for i in range(4):
        v1 = edges[i] / np.linalg.norm(edges[i])
        v2 = edges[(i + 1) % 4] / np.linalg.norm(edges[(i + 1) % 4])
        ang = np.degrees(np.arccos(np.clip(np.dot(v1, v2), -1, 1)))
        if abs(ang - 90) > angle_tol:
            return False

    return True


def contains_internal_atom(quad: Tuple[int, ...], terrace: List[int], positions: np.ndarray) -> bool:
    """Ensure fourfold ring candidate encloses no internal surface atoms."""
    coords2d = positions[list(quad)][:, :2]
    center = coords2d.mean(axis=0)
    angles = np.arctan2(coords2d[:, 1] - center[1], coords2d[:, 0] - center[0])
    order = np.argsort(angles)
    polygon = coords2d[order]

    path = MPLPath(polygon)
    quad_set = set(quad)

    for atom in terrace:
        if atom in quad_set:
            continue
        point = positions[atom][:2]
        if path.contains_point(point):
            return True

    return False


def prune_unique(candidates: List[Tuple[Tuple[int, ...], Tuple[int, ...]]]) -> List[Tuple[int, ...]]:
    """Symmetry prune site candidates to retain single equivalent representative."""
    unique_dict = {}
    for key, tri in candidates:
        if key not in unique_dict:
            unique_dict[key] = tri
    return list(unique_dict.values())


def get_layer_grouping(heights: np.ndarray, max_layer_span: float = 0.8) -> Tuple[np.ndarray, np.ndarray]:
    """Hierarchical 1D clustering on Z-heights to extract discrete atomic layers."""
    heights_2d = heights.reshape(-1, 1)
    Z = linkage(heights_2d, method='complete')
    clusters = fcluster(Z, t=max_layer_span, criterion='distance')
    
    unique_clusters = np.unique(clusters)
    layer_means = np.array([np.mean(heights[clusters == c]) for c in unique_clusters])
    
    # Sort cluster IDs by layer elevation ascending
    sorted_indices = np.argsort(layer_means)
    layer_means = layer_means[sorted_indices]
    
    # Remap cluster labels from bottom (1) to top (N)
    label_map = {old_label: new_label + 1 for new_label, old_label in enumerate(unique_clusters[sorted_indices])}
    remapped_clusters = np.array([label_map[c] for c in clusters])
    
    return remapped_clusters, layer_means


# Core Processing Engine
def analyze_slab(
    slab: Any,
    nn_dist: float,
    output_dir: Path,
    metadata_base: Dict[str, Any],
    site: list[int] | None = None,
    freeze_bottom: bool = False,
    freeze_fraction: float = 0.4
) -> Dict[str, Any]:
    """Enumerates surface sites directly on an ASE Atoms object."""
    positions = slab.get_positions()
    positions_frac = slab.get_scaled_positions()
    natoms = len(slab)
    heights = positions[:, 2]  # Direct Z-coordinate height assumption
    max_height = np.max(heights)

    # 1D Hierarchical Clustering for Layer Identification
    clusters, layer_means = get_layer_grouping(heights, max_layer_span=0.8)
    n_layers = len(layer_means)

    # Freeze bottom layers if requested or missing constraints
    if freeze_bottom and not slab.constraints:
        n_freeze = int(np.ceil(freeze_fraction * n_layers))
        frozen_clusters = set(range(1, n_freeze + 1))
        frozen_mask = np.isin(clusters, list(frozen_clusters))
        slab.set_constraint(FixAtoms(mask=frozen_mask))
    else:
        frozen_mask = np.array([False] * natoms)
        if slab.constraints:
            for constraint in slab.constraints:
                if isinstance(constraint, FixAtoms):
                    frozen_mask[constraint.index] = True
                    #frozen_mask |= constraint.index

    # Early exit for surface site override
    if site is not None:
        # Assemble metadata dictionary
        metadata = {
                **metadata_base,
                "detected_layer_count": int(n_layers),
                "cell": slab.get_cell().tolist(),
                "total_atoms": int(natoms),
                "surface_atoms": [int(i) for i in site],
                "frozen_atoms": [int(i) for i in np.where(frozen_mask)[0]]
                }

        output_dir.mkdir(parents=True, exist_ok=True)
        write(output_dir / "POSCAR", slab, vasp5=True, direct=True)
        (output_dir.parent / "surface_atoms.json").write_text(json.dumps(metadata, indent=2))

        return metadata

    # Neighbor calculations & undercoordination
    cutoff = 1.25 * nn_dist
    nl = NeighborList([cutoff] * natoms, self_interaction=False, bothways=True)
    nl.update(slab)

    coordination = np.array([len(nl.get_neighbors(i)[0]) for i in range(natoms)])
    bulk_coord = int(np.percentile(coordination, 90))
    surface_atoms_all = [i for i in range(natoms) if coordination[i] < bulk_coord and abs(heights[i] - max_height) < 2.5]

    # Symmetry setup
    sym_data = spglib.get_symmetry_dataset((slab.get_cell(), positions_frac, slab.get_atomic_numbers()), symprec=1e-3)
    equiv = sym_data.equivalent_atoms if sym_data is not None else list(range(natoms))

    unique_top = {}
    for i in surface_atoms_all:
        cls = equiv[i]
        if cls not in unique_top:
            unique_top[cls] = i
    pruned_surface_atoms = list(unique_top.values())

    # Bridge sites logic
    all_distances = [
        np.linalg.norm(positions[i] - positions[j])
        for i in range(natoms) for j in nl.get_neighbors(i)[0] if i < j
    ]
    bridge_cutoff = 1.2 * (np.min(all_distances) if all_distances else nn_dist)

    bridge_sites = []
    for i in surface_atoms_all:
        for j in nl.get_neighbors(i)[0]:
            if j in surface_atoms_all and i < j:
                if not crosses_pbc([i, j], positions_frac, threshold=0.5):
                    if np.linalg.norm(positions[i] - positions[j]) <= bridge_cutoff:
                        bridge_sites.append((i, j))

    bridge_sites_pruned = []
    seen_keys = set()
    for i, j in bridge_sites:
        key = tuple(sorted([equiv[i], equiv[j]]))
        if key not in seen_keys:
            seen_keys.add(key)
            bridge_sites_pruned.append((i, j))

    # Threefold and Fourfold logic
    max_edge = 1.15 * nn_dist

    # Select top layer atoms via cluster assignment
    top_cluster_id = n_layers
    top_layer_idx = np.where(clusters == top_cluster_id)[0]
    bulk_top_coord = int(np.median([coordination[i] for i in top_layer_idx])) if len(top_layer_idx) > 0 else 0
    valid_top_atoms = [i for i in top_layer_idx if coordination[i] == bulk_top_coord]

    # Group top atoms into terraces
    terrace_groups = [valid_top_atoms] if valid_top_atoms else []

    threefold_candidates = []
    fourfold_candidates = []

    for terrace in terrace_groups:
        terrace_coords = positions[terrace]
        if len(terrace) < 3:
            continue

        centroid = terrace_coords.mean(axis=0)
        cov = (terrace_coords - centroid).T @ (terrace_coords - centroid)
        eigvals, eigvecs = np.linalg.eigh(cov)
        local_normal = eigvecs[:, np.argmin(eigvals)]
        local_normal /= np.linalg.norm(local_normal)

        proj_coords = project_onto_plane(terrace_coords, local_normal)
        coords_2d = proj_coords[:, :2]

        try:
            tri = Delaunay(coords_2d)
            simplices = [tuple(terrace[v] for v in simplex) for simplex in tri.simplices]
            threefold_candidates.extend(simplices)
        except Exception:
            pass

        adjacency = build_terrace_adjacency(terrace, positions)
        for i in terrace:
            for j in adjacency.get(i, []):
                if j <= i:
                    continue
                for k in adjacency.get(j, []):
                    if k in (i, j) or k <= i:
                        continue
                    for l in adjacency.get(k, []):
                        if l in (i, j, k) or l <= i:
                            continue
                        if i in adjacency.get(l, []):
                            quad = tuple(sorted([i, j, k, l]))
                            if crosses_pbc(quad, positions_frac, threshold=0.8):
                                continue
                            if not is_rectangular_ring(quad, positions, heights, layer_tol=0.8):
                                continue
                            if contains_internal_atom(quad, terrace, positions):
                                continue
                            fourfold_candidates.append(quad)

    # Filter, classify, and cluster Threefolds
    filtered_threefolds = [
        tri_atoms for tri_atoms in threefold_candidates
        if is_connected_triangle(tri_atoms, nl)
        and not crosses_pbc(tri_atoms, positions_frac, threshold=0.8)
        and triangle_max_edge(tri_atoms, positions, max_edge)
    ]

    sub_layer_idx = np.where(clusters == (n_layers - 1))[0] if n_layers > 1 else np.array([], dtype=int)
    sub_positions = positions[sub_layer_idx]

    threefold_types = [classify_fcc_hcp(tri, positions, sub_positions) for tri in filtered_threefolds]
    threefold_sites = cluster_sites(filtered_threefolds, positions)
    threefold_types = [threefold_types[filtered_threefolds.index(t)] for t in threefold_sites]

    fcc_candidates = []
    hcp_candidates = []
    for tri, ttype in zip(threefold_sites, threefold_types):
        key = tuple(sorted([equiv[i] for i in tri]))
        if ttype == "fcc":
            fcc_candidates.append((key, tri))
        else:
            hcp_candidates.append((key, tri))

    fcc_sites_pruned = prune_unique(fcc_candidates)[:1]
    hcp_sites_pruned = prune_unique(hcp_candidates)[:1]

    # Filter and prune Fourfolds
    fourfold_candidates = list(set(fourfold_candidates))
    fourfold_unique = {}
    for quad in fourfold_candidates:
        key = tuple(sorted([equiv[i] for i in quad]))
        if key not in fourfold_unique:
            fourfold_unique[key] = quad

    fourfold_sites_pruned = list(fourfold_unique.values())[:1]

    # Assemble metadata dictionary
    metadata = {
        **metadata_base,
        "detected_layer_count": int(n_layers),
        "cell": slab.get_cell().tolist(),
        "total_atoms": int(natoms),
        "surface_atoms": [int(i) for i in pruned_surface_atoms],
        "bridge_sites": [[int(i), int(j)] for i, j in bridge_sites_pruned],
        "threefold_sites": [[int(i) for i in trip] for trip in hcp_sites_pruned + fcc_sites_pruned],
        "threefold_types": ["hcp"] * len(hcp_sites_pruned) + ["fcc"] * len(fcc_sites_pruned),
        "fourfold_sites": [[int(i) for i in quad] for quad in fourfold_sites_pruned],
        "fourfold_types": ["fourfold"] * len(fourfold_sites_pruned),
        "frozen_atoms": [int(i) for i in np.where(frozen_mask)[0]]
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    write(output_dir / "POSCAR", slab, vasp5=True, direct=True)
    (output_dir.parent / "surface_atoms.json").write_text(json.dumps(metadata, indent=2))

    return metadata


# High-Level Entry Points
def generate_slab_surface(
    element: str,
    miller: Tuple[int, int, int],
    output_dir: Path,
    packing: str = "fcc",
    lattconst: Optional[float] = None,
    target_thickness: float = 6.0,
    min_lateral_size: float = 8.0,
    vacuum: float = 15.0,
    freeze_fraction: float = 0.4
) -> Dict[str, Any]:
    """Generates POSCAR and metadata dictionary for a metallic surface slab."""
    bulk_met = bulk(element, packing, cubic=True) if lattconst is None else bulk(element, packing, a=lattconst, cubic=True)

    layers = 1
    while True:
        slab = surface(bulk_met, miller, layers=layers)
        slab.center(vacuum=vacuum, axis=2)
        normal = np.cross(slab.get_cell()[0], slab.get_cell()[1])
        normal /= np.linalg.norm(normal)
        heights = slab.get_positions() @ normal
        if (heights.max() - heights.min()) >= target_thickness:
            break
        layers += 1

    cell = slab.get_cell()
    rx = max(1, int(np.ceil(min_lateral_size / np.linalg.norm(cell[0]))))
    ry = max(1, int(np.ceil(min_lateral_size / np.linalg.norm(cell[1]))))
    slab = slab.repeat((rx, ry, 1))
    slab.set_pbc((True, True, True))

    bulk_dists = [bulk_met.get_distance(i, j, mic=True) for i in range(len(bulk_met)) for j in range(i + 1, len(bulk_met))]
    nn_dist = np.min(bulk_dists)

    metadata_base = {
        "element": element,
        "packing": packing,
        "miller": list(miller),
        "requested_layers": int(layers),
    }

    return analyze_slab(slab, nn_dist, output_dir, metadata_base, freeze_bottom=True, freeze_fraction=freeze_fraction)


def load_slab_from_poscar(
    poscar_path: Path,
    output_dir: Path,
    site: list[int] | None = None,
    freeze_bottom: bool = False,
    freeze_fraction: float = 0.4
) -> Dict[str, Any]:
    """Loads external POSCAR and performs surface site enumeration."""
    slab = read(poscar_path)

    # 1. Estimate nearest-neighbor distance safely across mixed-element systems
    positions = slab.get_positions()
    z_min = positions[:, 2].min()
    bulk_like_mask = positions[:, 2] < (z_min + 3.5)

    if np.sum(bulk_like_mask) > 1:
        dists = slab.get_all_distances(mic=True)[bulk_like_mask][:, bulk_like_mask]
        nn_dist = np.min(dists[dists > 0.5])
    else:
        # Fallback to global minimum pairwise distance
        dists = slab.get_all_distances(mic=True)
        nn_dist = np.min(dists[dists > 0.5])

    # 2. Preserve POSCAR elemental order (DO NOT use set())
    # dict.fromkeys preserves insertion order while removing duplicates
    ordered_elements = list(dict.fromkeys(slab.get_chemical_symbols()))

    metadata_base = {
        "source_file": str(poscar_path),
        "element": ordered_elements,
        "miller": None,
        "requested_layers": None,
    }

    return analyze_slab(
        slab, 
        nn_dist, 
        output_dir, 
        metadata_base, 
        site, 
        freeze_bottom=freeze_bottom, 
        freeze_fraction=freeze_fraction
    )

