function fp = fp_grey_calc_001(LATTICE, COORDINATES, numIons)

global ORG_STRUC
atomType = ORG_STRUC.atomType;

%fp_calc(LATTICE, COORDINATES, numIons) calculates fingerprint for
% structure defined by LATTICE, COORDINATES, numIons
numIons = sum(numIons); %Grey Fp
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, 1); %grey FP 
[none, fp, none] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons); %grey FP

end

