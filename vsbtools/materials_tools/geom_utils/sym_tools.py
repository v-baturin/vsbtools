from pymatgen.io.vasp import Poscar
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from collections import defaultdict

def extract_structure_info(poscar_path, element_order):
    """
    Extracts structural info from a POSCAR file.

    Args:
        poscar_path (str): Path to the POSCAR file.
        element_order (list of str): Desired order of elements in the output.

    Returns:
        dict: {
            "formula": str,
            "space_group": int,
            "inequivalent_positions": dict[str, int]
        }
    """
    structure = Poscar.from_file(poscar_path).structure
    sga = SpacegroupAnalyzer(structure, symprec=1e-3)
    symm_structure = sga.get_symmetrized_structure()

    # Count inequivalent positions
    inequivalent_counts = defaultdict(int)
    for group in symm_structure.equivalent_sites:
        element = group[0].specie.symbol
        inequivalent_counts[element] += 1

    # Format result with user-defined element order
    result = {
        "formula": structure.composition.reduced_formula,
        "space_group": sga.get_space_group_number(),
        "inequivalent_positions": {
            el: inequivalent_counts[el] for el in element_order if el in inequivalent_counts
        }
    }
    return result

if __name__ == "__main__":
    test_poscar = ("/home/vsbat/SYNC/00__WORK/2025-2026_MOLTEN_SALTS/STRUCTURE_GENERATION"
                   "/filtered_structures/B-Si-Mo/POSCAR_SiB2Mo5_agm003218039")
    print(extract_structure_info(test_poscar, ["Mo", "Si", "B"]))

