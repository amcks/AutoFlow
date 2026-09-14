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
@click.option('-s', '--slab', required=True, help='Slab element (e.g. Cu, Pt).')
@click.option('-m', '--miller', required=True, help='Comma-separated Miller indices (e.g. 1,1,1).')
@click.option('-a', '--adsorbate', required=True, help='Adsorbate SMILES string.')
@click.option('-l', '--lattconst', type=float, default=None, help='Lattice constant.')
@click.option('-p', '--packing', default='fcc', type=click.Choice(['fcc', 'hcp', 'bcc', 'bct']))
@click.option('-j', '--jobs', default=4, help='Max parallel screening jobs.')
@click.option('--generate-dft', is_flag=True, default=False, help='Prepare DFT single-point directories after screening.')
@click.option('--vasp-pp', type=click.Path(exists=True), default=None, help='Path to VASP POTCAR directory.')
@click.option('--mace-model', type=click.Path(exists=True), default=None, help='Path to MACE model checkpoint.')
def run_pipeline(slab, miller, adsorbate, lattconst, packing, jobs, generate_dft, vasp_pp, mace_model):
    """Execute standard screening pipeline."""
    # Defer heavy imports after click processing for lighter --help function
    from autoflow.pipeline import run_autoflow_pipeline
    from autoflow.generators.dft import generate_dft_input

    vasp_path = Path(vasp_pp) if vasp_pp else DEFAULT_VASP_POTENTIAL_PATH
    mace_path = Path(mace_model) if mace_model else DEFAULT_MACE_PATH
    miller_tuple = tuple(map(int, miller.split(',')))

    run_autoflow_pipeline(
        slab_element=slab,
        miller=miller_tuple,
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
    # Defer heavy imports after click processing for lighter --help function
    from autoflow.generators.dft import generate_dft_input

    generate_dft_input(screening_dir=Path(dir))


if __name__ == '__main__':
    main()
