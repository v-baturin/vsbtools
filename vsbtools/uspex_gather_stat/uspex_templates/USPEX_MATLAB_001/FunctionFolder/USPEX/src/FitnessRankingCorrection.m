function fitness = FitnessRankingCorrection(fitness)

global POP_STRUC
global ORG_STRUC

% ranking gives the order of fitness funtion values for this generation
[nothing, ranking] = sort(fitness);
POP_STRUC.ranking = ranking;
bad_rank = 0;
Step = length(ORG_STRUC.abinitioCode);

% all structures with distance less than toleranceFing are considered identical
% best structure is left, rest are moved to the end of the list

[ranking, bad_rank, fitness] = degradeSimilar(ORG_STRUC.toleranceFing, ranking, fitness, Step);

POP_STRUC.fitness = fitness;
POP_STRUC.ranking = ranking;
POP_STRUC.bad_rank = bad_rank;
safesave ('Current_POP.mat', POP_STRUC)
	
