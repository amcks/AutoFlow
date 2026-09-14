from pathlib import Path
from typing import List
from autoflow.core.config import DEFAULT_VASP_POTENTIAL_PATH

def write_kpoints(output_path: Path, comment: str, grid: str) -> None:
    """Writes a standard VASP KPOINTS file."""
    content = f"{comment}\n0\n{grid}\n0 0 0\n" if "Gamma" in grid or "Monkhorst" in grid else f"{comment}\n1\nreciprocal\n0.0  0.0  0.0  1.0\n"
    output_path.write_text(content)

def build_potcar(elements: List[str], output_path: Path, pot_base_dir: Path = DEFAULT_VASP_POTENTIAL_PATH) -> None:
    """Concatenates POTCAR files for a list of element symbols in order."""
    combined_content = bytearray()
    for el in elements:
        pot_file = pot_base_dir / el / "POTCAR"
        if not pot_file.exists():
            raise FileNotFoundError(f"POTCAR file for element '{el}' not found at {pot_file}")
        combined_content.extend(pot_file.read_bytes())
    
    output_path.write_bytes(combined_content)
