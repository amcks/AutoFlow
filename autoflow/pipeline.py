import sys
import subprocess
import shutil
from pathlib import Path
from ase.io import read

from autoflow.core.config import DEFAULT_VASP_POTENTIAL_PATH, DEFAULT_MACE_PATH, GAS_INCAR_TEMPLATE, SLAB_INCAR_TEMPLATE, SLURM_SCRIPT
from autoflow.core.utils import build_potcar, write_kpoints
from autoflow.generators.gas import generate_gas_phase
from autoflow.generators.inputs import generate_dockonsurf_input, generate_monoatomic_configurations
from autoflow.generators.slab import generate_slab_surface, load_slab_from_poscar
from autoflow.screening.post_analysis import run_post_analysis
from autoflow.screening.runner import run_parallel_screening


def run_autoflow_pipeline(
    smiles: str,
    max_parallel_jobs: int,
    slab_element: str | None = None,
    miller: tuple[int, int, int] | None = None,
    poscar_path: Path | None = None,
    latt_const: float | None = None,
    packing: str = "fcc",
    vasp_potential_path: Path = DEFAULT_VASP_POTENTIAL_PATH,
    mace_model_path: Path = DEFAULT_MACE_PATH,
) -> None:
    work_dir = Path.cwd()
    gas_dir = work_dir / "gas"
    slab_dir = work_dir / "slab"
    screening_dir = work_dir / "screening"

    # 1. Gas Phase Stage
    gas_elements = generate_gas_phase(smiles, gas_dir)
    build_potcar(gas_elements, gas_dir / "POTCAR")
    (gas_dir / "INCAR").write_text(GAS_INCAR_TEMPLATE)
    write_kpoints(gas_dir / "KPOINTS", "Gamma-point only", "Gamma")
    (gas_dir / "submit.sh").write_text(SLURM_SCRIPT)

    # 2. Slab Phase Stage
    slab_dir.mkdir(parents=True, exist_ok=True)
    
    if poscar_path is not None:
        # Defensive measure to resolve relative paths
        poscar_path = Path(poscar_path).expanduser().resolve()
        # Load user-provided POSCAR
        slab_meta = load_slab_from_poscar(
            poscar_path=poscar_path,
            output_dir=slab_dir,
        )
        # Ensure exact input POSCAR is copied into slab_dir
        shutil.copy(poscar_path, slab_dir / "POSCAR")
    else:
        # Generate slab via ASE
        slab_meta = generate_slab_surface(
            element=slab_element,
            miller=miller,
            output_dir=slab_dir,
            packing=packing,
            lattconst=latt_const,
        )

    # Handle single element string vs list of elements for alloys
    slab_elements = slab_meta["element"]
    if isinstance(slab_elements, str):
        slab_elements = [slab_elements]

    build_potcar(slab_elements, slab_dir / "POTCAR")
    (slab_dir / "INCAR").write_text(SLAB_INCAR_TEMPLATE)
    write_kpoints(slab_dir / "KPOINTS", "Slab k-points", "4 4 1")

    # 3. Enumeration Stage (inside screening/)
    combined_elements = slab_elements + gas_elements
    build_potcar(combined_elements, work_dir / "POTCAR")
    (work_dir / "INCAR").write_text(SLAB_INCAR_TEMPLATE)
    write_kpoints(work_dir / "KPOINTS", "Slab k-points", "4 4 1")
    (work_dir / "submit.sh").write_text(SLURM_SCRIPT)

    # Check adsorbate atom count
    n_atoms = len(read(gas_dir / "meta.xyz"))

    if n_atoms == 1:
        screening_dir.mkdir(exist_ok=True)
        generate_monoatomic_configurations(
            slab_path=slab_dir / "POSCAR",
            gas_path=gas_dir / "POSCAR",
            json_path=work_dir / "surface_atoms.json",
            output_dir=screening_dir,
        )
    else:
        inp_file = work_dir / "dockonsurf.inp"
        generate_dockonsurf_input(
            json_path=work_dir / "surface_atoms.json",
            gas_meta_path=gas_dir / "meta.xyz",
            output_file=inp_file,
        )

        dockonsurf_exec = shutil.which("dockonsurf.py")
        if not dockonsurf_exec:
            raise RuntimeError(
                "Executable 'dockonsurf.py' not found in system PATH. "
                "Please ensure DockonSurf is in your environment and PATH."
            )

        try:
            subprocess.run(
                ["dockonsurf.py", "-i", str(inp_file), "-f"],
                check=True,
                cwd=work_dir,
                capture_output=True,
                text=True
            )
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"DockonSurf failed with returncode {e.returncode}:\n{e.stderr}\n")
            raise e

    # File cleanup
    for_transport = ["dockonsurf.inp", "dockonsurf.log", "POSCAR", "POTCAR", "INCAR", "KPOINTS", "submit.sh", "surface_atoms.json"]
    screening_dir.mkdir(exist_ok=True)
    for file_name in for_transport:
        source = work_dir / file_name
        if source.is_file():
            shutil.move(source, screening_dir / file_name)

    # 4. Screening Stage
    run_parallel_screening(
        screening_dir=screening_dir,
        model_path=mace_model_path,
        max_jobs=max_parallel_jobs,
    )

    # Prerun cleanup
    for_removal = ["gfnff_topo", "gfnff_adjacency", "gfnff_lists.json"]
    for file_name in for_removal:
        source = work_dir / file_name
        source.unlink(missing_ok=True)

    # 5. Post-Analysis
    run_post_analysis(screening_dir)
