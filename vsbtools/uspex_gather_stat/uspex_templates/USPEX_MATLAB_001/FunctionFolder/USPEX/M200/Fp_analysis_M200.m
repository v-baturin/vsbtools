function Fp_analysis_M200(whichInd)
global ORG_STRUC
global POP_STRUC

LATTICE     = POP_STRUC.POPULATION(whichInd).LATTICE_2D;
numIons     = POP_STRUC.POPULATION(whichInd).numIons;
COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES_2D;
atomType    = ORG_STRUC.atomType;

[Ni, V, dist_matrix, typ_i, typ_j, ho, ht] = makeMatrices_2D(LATTICE, COORDINATES, numIons, atomType);
[order, FINGERPRINT, atom_fing] = fingerprint_calc_2D(Ni, V, dist_matrix, typ_i, typ_j, numIons, ho, ht);
POP_STRUC.POPULATION(whichInd).order =  order;
POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
POP_STRUC.POPULATION(whichInd).struc_entr = structureQuasiEntropy(numIons, atom_fing, numIons/sum(numIons));
POP_STRUC.POPULATION(whichInd).S_order    = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);
