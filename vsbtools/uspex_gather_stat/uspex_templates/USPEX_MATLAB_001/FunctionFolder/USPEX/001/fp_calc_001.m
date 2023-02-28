function fp = fp_calc_001(LATTICE, COORDINATES, numIons)

global ORG_STRUC
atomType = ORG_STRUC.atomType;

%fp_calc(LATTICE, COORDINATES, numIons) calculates fingerprint for
% structure defined by LATTICE, COORDINATES, numIons
% numIons = sum(numIons);%Grey init_Fp for varcomp
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType); %not grey FP for varcomp 001
[none, fp, none] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons); %not grey Fp For varcomp 001

end

