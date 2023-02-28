function [ranking, bad_rank, fitness] = META_CalcFitness(Step)

global POP_STRUC
global ORG_STRUC

fitness = [];
fitness = zeros(1,length(POP_STRUC.POPULATION));
for fit_loop = 1 : length(POP_STRUC.POPULATION)
  try
     fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(Step);
  catch
     fitness(fit_loop) = 1000000 + fit_loop;
  end
  if ORG_STRUC.varcomp == 1
     fitness(fit_loop) = fitness(fit_loop)/sum(POP_STRUC.POPULATION(fit_loop).numIons);
  end
end

[nothing, ranking] = sort(fitness);
if POP_STRUC.generation > 0
  [ranking, bad_rank, fitness] = degradeSimilar(ORG_STRUC.toleranceFing, ranking, fitness, Step);
else
   ranking = 1;
   bad_rank = 0;
end
POP_STRUC.ranking = ranking(1);   % to know which structure is the best!
