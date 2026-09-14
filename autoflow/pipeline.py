import sys
import subprocess
import shutil
from pathlib import Path
from ase.io import read

from autoflow.core.config import DEFAULT_VASP_POTENTIAL_PATH, DEFAULT_MACE_PATH, GAS_INCAR_TEMPLATE, SLAB_INCAR_TEMPLATE, SLURM_SCRIPT
from autoflow.core.utils import build_potcar, write_kpoints #,generate_slurm_script for later
from autoflow.generators.gas import generate_gas_phase
from autoflow.generators.inputs import generate_dockonsurf_input, generate_monoatomic_configurations
from autoflow.generators.slab import generate_slab_surface
from autoflow.screening.post_analysis import run_post_analysis
from autoflow.screening.runner import run_parallel_screening


def run_autoflow_pipeline(
    slab_element: str,
    miller: tuple[int, int, int],
    smiles: str,
    latt_const: float | None,
    packing: str,
    max_parallel_jobs: int,
    vasp_potential_path = DEFAULT_VASP_POTENTIAL_PATH,
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
    # Later change SLURM_SCRIPT to a function that takes in arguments like below
    #generate_slurm_script(gas_dir / "submit.sh", job_name="Gas", time="30:00:00")

    # 2. Slab Phase Stage
    slab_meta = generate_slab_surface(
        element=slab_element,
        miller=miller,
        output_dir=slab_dir,
        packing=packing,
        lattconst=latt_const,
    )
    # generate_slab_surface writes surface_atoms.json to slab_dir.parent (work_dir)
    build_potcar([slab_meta["element"]], slab_dir / "POTCAR")
    (slab_dir / "INCAR").write_text(SLAB_INCAR_TEMPLATE)
    write_kpoints(slab_dir / "KPOINTS", "Slab k-points", "4 4 1")

    # 3. Enumeration Stage (inside screening/)
    #screening_dir.mkdir(exist_ok=True)
    combined_elements = [slab_meta["element"]] + gas_elements
    build_potcar(combined_elements, work_dir / "POTCAR")
    (work_dir / "INCAR").write_text(SLAB_INCAR_TEMPLATE)
    write_kpoints(work_dir / "KPOINTS", "Slab k-points", "4 4 1")
    (work_dir / "submit.sh").write_text(SLURM_SCRIPT)
    # Later change SLURM_SCRIPT to a function that takes in arguments like below
    #generate_slurm_script(screening_dir / "submit.sh", job_name="Adsorb", time="12:00:00")

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

        # Execute DockonSurf using relative path references from work_dir
        dockonsurf_exec = shutil.which("dockonsurf.py")
        if not dockonsurf_exec:
            raise RuntimeError(
                    "Executable 'dockonsurf.py' not found in system PATH."
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
            #source.rename(screening_dir / file_name)

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
