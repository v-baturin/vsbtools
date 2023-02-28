function Fp_analysis_300(whichInd)
global ORG_STRUC
global POP_STRUC

LATTICE     = POP_STRUC.POPULATION(whichInd).LATTICE;
numIons     = POP_STRUC.POPULATION(whichInd).numIons;
COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
atomType    = ORG_STRUC.atomType;
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
[order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
POP_STRUC.POPULATION(whichInd).order       = order;
POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
POP_STRUC.POPULATION(whichInd).struc_entr  = structureQuasiEntropy(numIons, atom_fing, numIons/sum(numIons));
POP_STRUC.POPULATION(whichInd).S_order     = StructureOrder(FINGERPRINT, V, numIons, ...
                                             ORG_STRUC.deltaFing, ORG_STRUC.weight);
