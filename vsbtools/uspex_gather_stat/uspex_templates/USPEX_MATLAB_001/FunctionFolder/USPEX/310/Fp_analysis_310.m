function Fp_analysis_310(whichInd)

global ORG_STRUC
global POP_STRUC

LATTICE     = POP_STRUC.POPULATION(whichInd).LATTICE;
COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
numIons =     POP_STRUC.POPULATION(whichInd).numIons;
MtypeLIST =   POP_STRUC.POPULATION(whichInd).MtypeLIST;
numMols = sum(POP_STRUC.POPULATION(whichInd).numMols);
atomType = ORG_STRUC.atomType;

%This block is only used to calculate the Fp for Molcenters
coordinates = zeros(sum(numMols),3);
for i = 1 : sum(numMols)
    coordinates(i,:)=POP_STRUC.POPULATION(whichInd).MOLECULES(i).MOLCENTER/LATTICE;
end
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, coordinates, numMols, 1);
[order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numMols);
for i = 1: sum(numMols)
    POP_STRUC.POPULATION(whichInd).MOLECULES(i).order = order(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This block is only used to calculate the Fp of crystals, see Zhu, ActaB, 2012
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
Intra_map = Intra_MOL_dist(MtypeLIST, numIons, ORG_STRUC.STDMOL);
tmp = dist_matrix(1:sum(numIons), :);
tmp(find(Intra_map==0)) = 0;
dist_matrix(1:sum(numIons),:)  = tmp;

[order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);

POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
POP_STRUC.POPULATION(whichInd).struc_entr  = structureQuasiEntropy(numIons, atom_fing, numIons/sum(numIons));
POP_STRUC.POPULATION(whichInd).S_order     = StructureOrder(FINGERPRINT, V, numIons, ...
                                             ORG_STRUC.deltaFing, ORG_STRUC.weight);
POP_STRUC.POPULATION(whichInd).order       = order;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
