function Fp_analysis_301(whichInd)
global ORG_STRUC
global POP_STRUC

LATTICE = POP_STRUC.POPULATION(whichInd).LATTICE;
numIons = sum(POP_STRUC.POPULATION(whichInd).numIons); %!Grey Fp for varcomp
COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, 1); 
[order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons); 
POP_STRUC.POPULATION(whichInd).order       =  order;
POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
POP_STRUC.POPULATION(whichInd).struc_entr  = structureQuasiEntropy(numIons, atom_fing, 1);
POP_STRUC.POPULATION(whichInd).S_order     = StructureOrder(FINGERPRINT, V, numIons, ...
                                             ORG_STRUC.deltaFing, ORG_STRUC.weight);
