from functools import lru_cache
from ..crystal_entry import CrystalEntry
from USPEX.components import Atomistic
from USPEX.DataModel.Engine import Engine
from USPEX.DataModel.Flavour import Flavour
from USPEX.DataModel.Entry import Entry
import numpy as np
from USPEX.Atomistic.RadialDistributionUtility import RadialDistributionUtility, TOLERANCE_DEFAULT

Engine.createEngine(":memory:")
atomistic = Atomistic()


class USPEXBridge:
    def __init__(self, elements, legacy, tol_FP=None):
        tol_FP = tol_FP or TOLERANCE_DEFAULT
        self.rdu = RadialDistributionUtility(symbols=elements, suffix='origin', legacy=legacy, tolerance=tol_FP)
        self.uspex_entry_extensions = dict(atomistic=(atomistic, atomistic.propertyExtension.propertyTable),
                          radialDistributionUtility=(self.rdu, self.rdu.propertyExtension.propertyTable))
        self.id=-1

    @lru_cache(maxsize=1000)
    def uspex_entry_from_de(self, de_entry: CrystalEntry) -> "Entry":
        types, coords, cell = ([s.species_string for s in de_entry.structure],
                               de_entry.structure.cart_coords,
                               de_entry.structure.lattice.matrix)
        uspex_structure = atomistic.AtomicStructureRepresentation.structureType(atomTypes=types, coordinates=coords,
                                                                                    cell=cell)
        self.id += 1
        return Entry.newEntry(Flavour(extensions=self.uspex_entry_extensions,
                                      **{'.howCome': 'Seeds', '.parent': None, '.label': self.id},
                                      **atomistic.atomicDisassemblerType(
                                          np.arange(len(uspex_structure)).reshape((-1, 1))).disassemble(
                                          uspex_structure)))

    def fp_dist(self, de_entry_1: CrystalEntry, de_entry_2: CrystalEntry) -> float:
        return self.rdu.dist(self.uspex_entry_from_de(de_entry_1), self.uspex_entry_from_de(de_entry_2))
