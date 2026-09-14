import json
import shutil
from pathlib import Path
from typing import Dict, Any, List, Optional
from ase.io import read, write

from autoflow.core.config import SLAB_INCAR_TEMPLATE, SLURM_SCRIPT
from autoflow.core.utils import write_kpoints


def generate_dft_input(
    screening_dir: Path,
    output_dir: Optional[Path] = None,
    template_dir: Optional[Path] = None,
) -> List[Path]:
    """Reads dft_selection.json and prepares DFT calculation subdirectories.

    Parameters
    ----------
    screening_dir : Path
        Directory containing dft_selection.json and conf_* subdirectories.
    output_dir : Optional[Path]
        Target directory to output DFT calculation folders (defaults to screening_dir / "DFT").
    template_dir : Optional[Path]
        Directory containing custom INCAR, KPOINTS, POTCAR, or submit.sh templates.

    Returns
    -------
    List[Path]
        List of generated DFT run directory paths.
    """
    dft_json = screening_dir / "dft_selection.json"
    if not dft_json.exists():
        raise FileNotFoundError(f"Selection file not found at {dft_json}")

    dft_data = json.loads(dft_json.read_text())
    selection: List[Dict[str, Any]] = dft_data.get("selection", [])

    if output_dir is None:
        output_dir = screening_dir / "DFT"
    output_dir.mkdir(parents=True, exist_ok=True)

    prepared_dirs: List[Path] = []

    for idx, entry in enumerate(selection):
        conf_name = entry["conf"]
        method = entry["method"]

        conf_dir = screening_dir / conf_name
        xyz_file = conf_dir / f"relaxed_{method}.xyz"

        if not xyz_file.exists():
            raise FileNotFoundError(f"Structure file not found: {xyz_file}")

        atoms = read(xyz_file)
        run_dir = output_dir / f"{idx}"
        #run_dir = output_dir / f"run_{idx}"
        run_dir.mkdir(parents=True, exist_ok=True)

        # Write POSCAR
        write(run_dir / "POSCAR", atoms, format="vasp")

        # Handle Template/VASP Support Files
        for fname in ["INCAR", "KPOINTS", "POTCAR", "submit.sh"]:
            src = template_dir / fname if template_dir else screening_dir / fname
            dst = run_dir / fname

            if src.exists():
                shutil.copy(src, dst)
            else:
                # Fallback to defaults defined in core configs
                if fname == "INCAR":
                    dst.write_text(SLAB_INCAR_TEMPLATE)
                elif fname == "KPOINTS":
                    write_kpoints(dst, "DFT K-Points", "4 4 1")
                elif fname == "submit.sh":
                    dst.write_text(SLURM_SCRIPT)
                elif fname == "POTCAR" and (screening_dir / "POTCAR").exists():
                    shutil.copy(screening_dir / "POTCAR", dst)
                else:
                    dst.write_text(f"# Placeholder for {fname}\n")

        prepared_dirs.append(run_dir)

    return prepared_dirs
