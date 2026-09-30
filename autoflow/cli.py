import click
from pathlib import Path
from autoflow.core.config import (
    DEFAULT_VASP_POTENTIAL_PATH,
    DEFAULT_MACE_PATH,
)
from autoflow.core.utils import (
        ATOM_INDEX_PARAM,
        MILLER_INDEX_PARAM,
)

@click.group()
def main():
    """Autoflow: Automated surface adsorption screening pipeline."""
    pass


@main.command(name="run")
@click.option('-s', '--slab', default=None, help='Slab element (e.g. Cu, Pt). Required unless --poscar is used.')
@click.option('-m', '--miller', type=MILLER_INDEX_PARAM, default=None, help='Comma-separated Miller indices (e.g. 1,1,1). Required unless --poscar is used.')
@click.option('--poscar-slab', type=click.Path(exists=True, path_type=Path), default=None, help='Path to optional slab POSCAR file.')
@click.option('--poscar-gas', type=click.Path(exists=True, path_type=Path), default=None, help='Path to optional gas POSCAR file.')
@click.option('--site-slab', type=ATOM_INDEX_PARAM, default=None, help='Comma-separated atomic indices or grouped tuples (e.g. "43,5,(10,22)") for surface site override in supplied slab POSCAR file.')
@click.option('--site-gas', type=ATOM_INDEX_PARAM, default=None, help='Comma-separated atomic indices or grouped tuples for adsorbate molecule anchor points override in supplied gas POSCAR file.')
@click.option('-a', '--adsorbate', default=None, help='Adsorbate SMILES string.')
@click.option('-l', '--lattconst', type=float, default=None, help='Lattice constant.')
@click.option('-p', '--packing', default='fcc', type=click.Choice(['fcc', 'hcp', 'bcc', 'bct']))
@click.option('-j', '--jobs', default=4, help='Max parallel screening jobs.')
@click.option('--vasp-pp', type=click.Path(exists=True), default=None, help='Path to VASP POTCAR directory if not specified via bash variable.')
@click.option('--mace-model', type=click.Path(exists=True), default=None, help='Path to MACE model file if not specified via bash variable.')
@click.option('--no-screen', is_flag=True, default=False, help='[Flag] Option to stop after adsorption mode enumeration without performing MLIP screening.')
@click.option('--generate-dft', is_flag=True, default=False, help='[Flag] Option to prepare DFT single-point directories after screening.')
def run_pipeline(slab, miller, poscar_slab, poscar_gas, site_slab, site_gas, adsorbate, lattconst, packing, jobs, vasp_pp, mace_model, no_screen, generate_dft):
    """Execute standard screening pipeline."""
    # Enforce slab input mutual exclusivity
    if poscar_slab is None:
        if not slab or not miller:
            raise click.UsageError("Must provide both --slab and --miller, or provide an existing --poscar-slab file.")
        if site_slab is not None:
            raise click.UsageError("Cannot use --site-slab without providing --poscar-slab.")
    elif slab or miller:
        raise click.UsageError("Cannot specify --slab or --miller when using --poscar-slab.")

    # Enforce gas input mutual exclusivity
    if poscar_gas is None:
        if not adsorbate:
            raise click.UsageError("Must provide --adsorbate SMILES, or provide an existing --poscar-gas file.")
    elif adsorbate:
        raise click.UsageError("Cannot specify --adsorbate when using --poscar-gas.")
    elif site_gas is None:
        raise click.UsageError("Must specify --site-gas when using --poscar-gas.")

    from autoflow.pipeline import run_autoflow_pipeline
    from autoflow.generators.dft import generate_dft_input

    vasp_path = Path(vasp_pp) if vasp_pp else DEFAULT_VASP_POTENTIAL_PATH
    mace_path = Path(mace_model) if mace_model else DEFAULT_MACE_PATH

    run_autoflow_pipeline(
        slab_element=slab,
        miller=miller,
        poscar_slab=poscar_slab,
        poscar_gas=poscar_gas,
        site_slab=site_slab,
        site_gas=site_gas,
        smiles=adsorbate,
        latt_const=lattconst,
        packing=packing,
        max_parallel_jobs=jobs,
        vasp_potential_path=vasp_path,
        mace_model_path=mace_path,
        run_screening=not no_screen,
    )

    if generate_dft:
        if no_screen:
            raise click.UsageError("Cannot use --generate-dft alongside --no-screen as screening is required for DFT inputs generation.")
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
