### FIRST_STEPS - SIMPLE RELAXATION WITH MATTERGEN/MATTERSIM VENV
from typing import Tuple
from pathlib import Path
import numpy as np
from ase import Atoms
from ase.io import read
from mattersim.forcefield.potential import MatterSimCalculator
from mattersim.applications.relax import Relaxer

sio2_flavor_paths = dict(
mg=Path("2442POSCAR"),
agm=Path("agm002228342POSCAR"),
oqmd=Path("oqmd_104526POSCAR")
)

flavor = "mg"


# load the structure
flavor_path = sio2_flavor_paths[flavor]
ats: Atoms = read(flavor_path)

def relax(ats):


    # Distort the positions
    ats.positions += 0.1 * np.random.randn(len(ats), 3)

    #attach the calculator to the Atoms object
    ats.calc = MatterSimCalculator()

    # Initialize the relaxation object
    relaxer = Relaxer(
        optimizer='BFGS',
        filter='ExpCellFilter',
        constrain_symmetry=False
    )

    # Run the relaxation
    flag, result = relaxer.relax(ats, steps=500)

    print(f"{'relaxed' if flag else 'relaxation failed'}")
    return result
# result.write(f"{flavor_path.name}_relaxed", format='vasp', vasp5=True)



# import unittest
# from pathlib import Path
# from ase.io import read
# from ..coordination_tools import compute_average_coordination_atoms, _to_torch_from_atoms, _soft_neighbor_counts_per_A, compute_target_share_atoms
#
# PATH_WITH_TESTS = Path(__file__).parent
#
# class symmetry_tools_Test(unittest.TestCase):
#     def setUp(self):
#         self.system = read(PATH_WITH_TESTS / "Na3B7O12_POSCAR")
#
#     def test_pretty_phasediag(self):
#         self.assertAlmostEqual(compute_average_coordination_atoms(self.system, type_A="B", type_B="O", kernel='sigmoid'),
#                                3.43, delta=0.1)
#         print(_soft_neighbor_counts_per_A(*_to_torch_from_atoms(self.system), type_A="B", type_B="O"))
#         self.assertAlmostEqual(compute_target_share_atoms(self.system, type_A="B", type_B="O", target=3.),
#                                0.57, delta=0.05)
