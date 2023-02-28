function Fp_analysis_001(whichInd)

global ORG_STRUC
global POP_STRUC

LATTICE     = POP_STRUC.POPULATION(whichInd).LATTICE;
numIons     = POP_STRUC.POPULATION(whichInd).numIons; %!not grey Fp for varcomp 001
COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
atomType = ORG_STRUC.atomType;

[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType); %not grey Fp for varcomp

[order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons); %not grey Fp For varcomp

POP_STRUC.POPULATION(whichInd).order       =  order;
POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
POP_STRUC.POPULATION(whichInd).struc_entr  = structureQuasiEntropy(numIons, atom_fing, numIons/sum(numIons));
weight = CalcWeight_001(numIons);
POP_STRUC.POPULATION(whichInd).S_order     = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, weight);


