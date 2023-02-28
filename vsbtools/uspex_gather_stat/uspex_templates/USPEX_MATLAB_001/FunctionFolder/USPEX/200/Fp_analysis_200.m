function Fp_analysis_200(whichInd)
global ORG_STRUC
global POP_STRUC

LATTICE = POP_STRUC.POPULATION(whichInd).Surface_LATTICE;
COORDINATES = POP_STRUC.POPULATION(whichInd).Surface_COORDINATES;
numIons = POP_STRUC.POPULATION(whichInd).Surface_numIons;
[Ni, V, dist_matrix, typ_i, typ_j, ho, ht] = makeMatrices_2D...
               (LATTICE, COORDINATES, numIons, ORG_STRUC.atomType);
[order, FINGERPRINT, atom_fing] = fingerprint_calc_2D...
               (Ni, V, dist_matrix, typ_i, typ_j, numIons, ho, ht);
POP_STRUC.POPULATION(whichInd).Surface_order =  order;
POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
POP_STRUC.POPULATION(whichInd).struc_entr = structureQuasiEntropy(numIons, atom_fing, 1);
POP_STRUC.POPULATION(whichInd).S_order    = StructureOrder...
     (FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);

