import ast
import click
from pathlib import Path
from typing import List
from autoflow.core.config import DEFAULT_VASP_POTENTIAL_PATH

class MillerIndexParam(click.ParamType):
    """
    Click parameter type that converts strings like '1,1,1' or '1, 0, 0'
    strictly into a 3-element tuple of integers like (1,1,1).
    """
    name = "miller_indices"

    def convert(self, value, param, ctx):
        if value is None:
            return None

        if isinstance(value, (tuple, list)) and len(value) == 3 and all(isinstance(x, int) for x in value):
            return tuple(value)

        cleaned = str(value).strip()
        if not cleaned:
            return None

        try:
            parsed = [int(x.strip()) for x in cleaned.split(",") if x.strip()]
        except ValueError:
            self.fail(
                f"'{value}' is not a valid list of Miller indices. "
                f"Expected 3 comma-separated integers, e.g., '1,1,1'.",
                param,
                ctx,
            )

        if len(parsed) != 3:
            self.fail(
                f"Miller indices must consist of exactly 3 integers, got {len(parsed)} in '{value}'.",
                param,
                ctx,
            )

        return tuple(cleaned)

class AtomIndexParam(click.ParamType):
    """
    Click parameter type that converts strings like '1,4,(5,6)' or '1, 4, (5,6,7)'
    into Python structures: [1, 4, (5, 6)] or [1, 4, (5, 6, 7)].
    """
    name = "atom_indices"

    def convert(self, value, param, ctx):
        if value is None:
            return None

        # Pass through if it's already a list/tuple (e.g. called programmatically)
        if isinstance(value, (list, tuple)):
            return value

        cleaned = str(value).strip()
        if not cleaned:
            return []

        # Wrap in brackets if not already enclosed, turning "1,4,(5,6)" into "[1,4,(5,6)]"
        if not (cleaned.startswith("[") and cleaned.endswith("]")):
            expr = f"[{cleaned}]"
        else:
            expr = cleaned

        try:
            parsed = ast.literal_eval(expr)
        except (ValueError, SyntaxError):
            self.fail(
                f"'{value}' is not a valid index specification. "
                f"Expected integers or tuple groupings, e.g., '1,4,(5,6)'.",
                param,
                ctx,
            )

        if not isinstance(parsed, (list, tuple)):
            self.fail(f"Invalid format: '{value}'", param, ctx)

        # Validate elements inside the container
        validated = []
        for elem in parsed:
            if isinstance(elem, int):
                validated.append(elem)
            elif isinstance(elem, tuple) and all(isinstance(x, int) for x in elem):
                validated.append(elem)
            else:
                self.fail(
                    f"Invalid element '{elem}' in '{value}'. "
                    f"Elements must be integers or tuples of integers.",
                    param,
                    ctx,
                )

        return list(validated)

# Reusable instance for cli.py
MILLER_INDEX_PARAM = MillerIndexParam()
ATOM_INDEX_PARAM = AtomIndexParam()

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

