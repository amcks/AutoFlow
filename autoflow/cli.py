import click
from pathlib import Path
from autoflow.core.config import (
    DEFAULT_VASP_POTENTIAL_PATH,
    DEFAULT_MACE_PATH,
)

@click.group()
def main():
    """Autoflow: Automated surface adsorption screening pipeline."""
    pass


@main.command(name="run")
@click.option('-s', '--slab', default=None, help='Slab element (e.g. Cu, Pt). Required unless --poscar is used.')
@click.option('-m', '--miller', default=None, help='Comma-separated Miller indices (e.g. 1,1,1). Required unless --poscar is used.')
@click.option('--poscar', type=click.Path(exists=True, path_type=Path), default=None, help='Path to optional user-supplied POSCAR file.')
@click.option('-a', '--adsorbate', required=True, help='Adsorbate SMILES string.')
@click.option('-l', '--lattconst', type=float, default=None, help='Lattice constant.')
@click.option('-p', '--packing', default='fcc', type=click.Choice(['fcc', 'hcp', 'bcc', 'bct']))
@click.option('-j', '--jobs', default=4, help='Max parallel screening jobs.')
@click.option('--generate-dft', is_flag=True, default=False, help='Prepare DFT single-point directories after screening.')
@click.option('--vasp-pp', type=click.Path(exists=True), default=None, help='Path to VASP POTCAR directory.')
@click.option('--mace-model', type=click.Path(exists=True), default=None, help='Path to MACE model checkpoint.')
def run_pipeline(slab, miller, poscar, adsorbate, lattconst, packing, jobs, generate_dft, vasp_pp, mace_model):
    """Execute standard screening pipeline."""
    # Enforce input mutual exclusivity
    if poscar is None:
        if not slab or not miller:
            raise click.UsageError("You must provide both --slab and --miller, or provide an existing --poscar file.")
    elif slab or miller:
        raise click.UsageError("Cannot specify --slab or --miller when using --poscar.")

    from autoflow.pipeline import run_autoflow_pipeline
    from autoflow.generators.dft import generate_dft_input

    vasp_path = Path(vasp_pp) if vasp_pp else DEFAULT_VASP_POTENTIAL_PATH
    mace_path = Path(mace_model) if mace_model else DEFAULT_MACE_PATH
    miller_tuple = tuple(map(int, miller.split(','))) if miller else None

    run_autoflow_pipeline(
        slab_element=slab,
        miller=miller_tuple,
        poscar_path=poscar,
        smiles=adsorbate,
        latt_const=lattconst,
        packing=packing,
        max_parallel_jobs=jobs,
        vasp_potential_path=vasp_path,
        mace_model_path=mace_path,
    )

    if generate_dft:
        screening_dir = Path.cwd() / "screening"
        generate_dft_input(screening_dir=screening_dir)


@main.command(name="prep-dft")
@click.option('-d', '--dir', type=click.Path(exists=True), default="screening", help="Screening directory path.")
def prep_dft(dir):
    """Standalone command to generate DFT inputs."""
    from autoflow.generators.dft import generate_dft_input

    generate_dft_input(screening_dir=Path(dir))


if __name__ == '__main__':
    main()
