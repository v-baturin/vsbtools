function fitness = CalcFitness_M300()

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));

for i = 1:length(POP_STRUC.POPULATION)
     factor = POP_STRUC.POPULATION(i).numIons/ORG_STRUC.numIons;    
     fitness(i) = POP_STRUC.POPULATION(i).Enthalpies(end)/factor;
end

