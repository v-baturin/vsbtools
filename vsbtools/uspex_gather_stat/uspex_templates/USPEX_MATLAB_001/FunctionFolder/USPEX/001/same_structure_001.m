function same = same_structure_001(ind_pop, ind_cl)

global ORG_STRUC
global POP_STRUC
global CLUSTERS

same = 1;

numIons1 = POP_STRUC.POPULATION(ind_pop).numIons;
numIons2 = CLUSTERS.Structure(ind_cl).numIons;

%comparison of compositions
if ~isequal(numIons1,numIons2)
    same = 0;
end

%comparison of energies
toleranceE = 0.15;
if same == 1
    enthalpy1 = POP_STRUC.POPULATION(ind_pop).Enthalpies(end);
    enthalpy2 = CLUSTERS.Structure(ind_cl).enthalpy;
    if abs(enthalpy1-enthalpy2) > toleranceE
        same = 0;
    end
end

%comparison of fingerprints
toleranceF = ORG_STRUC.toleranceFing;
if same == 1
   fp1 = POP_STRUC.POPULATION(ind_pop).FINGERPRINT;
   fp2 = CLUSTERS.Structure(ind_cl).fingerprint;
   weight = CalcWeight_001(numIons1);
   dist_fp = cosineDistance(fp1, fp2, weight);
   if dist_fp > toleranceF
       same = 0;
   end
end