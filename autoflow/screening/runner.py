import json
import gc
import logging
from pathlib import Path
from typing import Optional, Dict, Any, List
import numpy as np
from ase.io import read, write
from ase.constraints import FixAtoms
from ase.optimize import FIRE, LBFGS
from xtb.ase.calculator import XTB
from mace.calculators import MACECalculator
from concurrent.futures import ProcessPoolExecutor, as_completed

logger = logging.getLogger(__name__)

def run_single_relaxation(conf_dir: Path, model_path: Path) -> Optional[Dict[str, Any]]:
    """Runs GFN-FF prerelaxation followed by MACE-FT relaxation for a single configuration directory."""
    poscar_path = conf_dir / "POSCAR"
    surface_json = conf_dir.parent / "surface_atoms.json"
    
    if not poscar_path.exists() or not surface_json.exists():
        print("here")
        return None

    atoms = read(poscar_path)
    slab_info = json.loads(surface_json.read_text())
    n_slab = slab_info["total_atoms"]

    # Freeze slab atoms
    mask = np.ones(len(atoms), dtype=bool)
    mask[:n_slab] = False
    atoms.set_constraint(FixAtoms(mask=~mask))

    try:
        # Pre-relaxation with GFN-FF
        atoms.set_pbc(False)
        atoms.calc = XTB(method="GFN-FF", charge=0, spin=0, maxiter=250, electronic_temperature=3000)
        opt_fire = FIRE(atoms, logfile=str(conf_dir / "prerun.log"))
        opt_fire.run(fmax=2, steps=10)
        
        atoms_relaxed = atoms.copy()
        atoms_relaxed.set_pbc(True)

        # Main Ensemble Relaxation (MACE-FT)
        mace_calc = MACECalculator(model_paths=str(model_path), device="cpu", default_dtype="float64")
        atoms_relaxed.calc = mace_calc
        
        opt_lbfgs = LBFGS(atoms_relaxed, logfile=str(conf_dir / "MACE-FT.log"))
        opt_lbfgs.run(fmax=0.2, steps=100)

        energy = float(atoms_relaxed.get_potential_energy())
        max_force = float(np.linalg.norm(atoms_relaxed.get_forces(), axis=1).max())

        atoms_relaxed.info["method"] = "MACE-FT"
        atoms_relaxed.info["energy"] = energy
        atoms_relaxed.info["max_force"] = max_force
        atoms_relaxed.calc = None

        write(conf_dir / "relaxed_MACE-FT.xyz", atoms_relaxed)
        
        res = [{"method": "MACE-FT", "energy": energy, "max_force": max_force}]
        (conf_dir / "ensemble_screen.json").write_text(json.dumps(res, indent=2))
        
        del mace_calc
        gc.collect()
        return res

    except Exception as e:
        (conf_dir / "FAILED").write_text(str(e))
        return None

def _worker_wrapper(conf_dir: Path, model_path: Path) -> tuple[Path, Optional[List[Dict[str, Any]]]]:
    """Helper wrapper to pair configuration path with execution result."""
    res = run_single_relaxation(conf_dir, model_path)
    return conf_dir, res

def run_parallel_screening(
    screening_dir: Path,
    model_path: Path,
    max_jobs: int = 4
) -> Dict[str, Optional[List[Dict[str, Any]]]]:
    """
    Executes MLIP relaxation across all conf_* subdirectories in parallel.
    
    Parameters
    ----------
    screening_dir : Path
        Directory containing configuration subdirectories (conf_0, conf_1, ...).
    model_path : Path
        Path to the trained interatomic potential checkpoint.
    max_jobs : int
        Maximum number of concurrent process workers.
    """
    conf_dirs = sorted([d for d in screening_dir.glob("conf_*") if d.is_dir()])
    
    if not conf_dirs:
        logger.warning(f"No configuration directories found in {screening_dir}")
        return {}

    results: Dict[str, Optional[List[Dict[str, Any]]]] = {}

    # Force sequential execution if max_jobs == 1 (easier debugging)
    if max_jobs == 1:
        logger.info("Executing screening sequentially (max_jobs=1)...")
        for d in conf_dirs:
            _, res = _worker_wrapper(d, model_path)
            results[d.name] = res
        return results

    logger.info(f"Starting parallel relaxation for {len(conf_dirs)} configurations across {max_jobs} workers...")

    # Parallel execution via ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=max_jobs) as executor:
        # Submit tasks to executor pool
        future_to_conf = {
            executor.submit(_worker_wrapper, conf_dir, model_path): conf_dir
            for conf_dir in conf_dirs
        }

        for future in as_completed(future_to_conf):
            conf_dir = future_to_conf[future]
            try:
                _, res = future.result()
                results[conf_dir.name] = res
                if res is not None:
                    logger.info(f"Completed relaxation: {conf_dir.name}")
                else:
                    logger.error(f"Relaxation failed or skipped: {conf_dir.name}")
            except Exception as exc:
                logger.error(f"Unhandled exception in process for {conf_dir.name}: {exc}")
                results[conf_dir.name] = None

    return results
