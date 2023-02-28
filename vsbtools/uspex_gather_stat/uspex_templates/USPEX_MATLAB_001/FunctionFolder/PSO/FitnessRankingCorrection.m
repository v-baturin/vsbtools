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

% for purely softmutated generations we do a clusterization starting from second generation
% best 90% of the structures are clustered into approx. populationSize clusters
% only best structure per cluster is taken.
% For this we tune the tolerance to have approx. populationSize structures before badrank
% Tuning is done via binary search algorithm

POP_STRUC.fitness = fitness;
POP_STRUC.ranking = ranking;
POP_STRUC.bad_rank = bad_rank;
save ('Current_POP.mat', 'POP_STRUC')

