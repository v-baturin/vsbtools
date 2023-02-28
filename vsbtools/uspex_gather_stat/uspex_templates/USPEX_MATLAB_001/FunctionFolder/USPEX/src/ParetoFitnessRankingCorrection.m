function fitness = ParetoFitnessRankingCorrection(fitness, Enthalpies)
global POP_STRUC
global ORG_STRUC

paretoRanking = ORG_STRUC.paretoRanking;
LENFIT = length(fitness);
Properties = [];
ErroR = 0 ;
originalFitness = fitness;
%% first we rank properties with respect to one property and then we take the structures which 
%% pass limitations and rank them according to two properties, and then put the structures 
%% which couldn't pass the the limitations, at the end of this structures in ranking.

%% First normal ranking according to fitness
[nothing, firstRank]  = sort(fitness);
firstRank = firstRank';
ParetoF1(:,1) = 1 : LENFIT;
ParetoFront = ParetoF1;

bad_rank = 0;
Step = length(ORG_STRUC.abinitioCode);

% all structures with distance less than toleranceFing are considered identical
% best structure is left, rest are moved to the end of the list
[Ranking, bad_rank, fitness] = degradeSimilar(ORG_STRUC.toleranceFing, firstRank, fitness, Step);

%%% to keep the repetitive structures at the end
for a = 1 : length(Ranking)
    if Ranking(a) == firstRank(end)
	sameStructures = Ranking(a+1:end);
	break;
    end
end
for a = 1 : length(sameStructures)
    fitness(sameStructures(a)) = 100000;
end
%%-----------------------------
firstStep = [Ranking, ParetoFront];


if paretoRanking ~= 1    %% paretoRanking > 1  : having more than 2 property to optimize simultaneously
    Properties = CalcPropforPareto(paretoRanking);
end
ii = 0;
Fit = []; ENE = []; PROP = []; ID = [];
for i = 1 : LENFIT
    if Enthalpies(i) < 1 && abs(fitness(i)) < 9999  
        ii = ii + 1 ;
        ID(ii)   = i;
        ENE(ii)  = Enthalpies(i);
        Fit(ii)  = fitness(i);
        if paretoRanking ~= 1
  	    PROP(ii) = Properties(i);
        end	
    end
end
%% Second Rank according to 2 or 3 properties
if paretoRanking == 1
    [secondRank, Dist, ParetoF2] = ParetoRanking_2D(Fit, ENE, ID);    % IDs are the number of Good Structures
    secondStep = [secondRank, ParetoF2];
    [Ranking, ParetoFront] = RANKING(secondStep, firstStep);
else
    [secondRank, Dist, ParetoF2] = ParetoRanking_3D(Fit, ENE, PROP, ID);
    secondStep = [secondRank, ParetoF2];
    [Ranking, ParetoFront] = RANKING(secondStep, firstStep);
end
initialStep = [Ranking, ParetoFront];
%% Third: Only for Dielectric constant, again ranking 
%% with respect to 2 or 3 properties, this time we have one 
%% more constrain, if no structures satisfy these constrains
%% ranking will not be performed, and previous rank is used.
if ORG_STRUC.optType == 6
    Fit = []; ENE = []; ID = []; PROP = [];
    ii      = 0;
    for i = 1 : LENFIT
       if Enthalpies(i) < 1 && abs(fitness(i)) < 9999 && fitness(i) > 1 && POP_STRUC.POPULATION(i).hardness > 2   %% Here we don't consider stability
            ii = ii + 1 ;
            ID(ii)   = i;
            ENE(ii)  = Enthalpies(i);
            Fit(ii)  = fitness(i);
            if paretoRanking ~= 1
                PROP(ii) = Properties(i);
            end
        end
    end
    if ii ~= 0
        if paretoRanking == 1
	    [thirdRank, Dist, ParetoF3] = ParetoRanking_2D(Fit, ENE, ID);
            thirdStep = [thirdRank, ParetoF3];
            [Ranking, ParetoFront] = RANKING(thirdStep, initialStep);
        else
            [thirdRank, Dist, ParetoF3] = ParetoRanking_3D(Fit, ENE, PROP, ID);
            thirdStep = [thirdRank, ParetoF3];
            [Ranking, ParetoFront] = RANKING(thirdStep, initialStep);
        end
    else
	disp('This system cannot satisfy hardness > 2GPa condition, ranking is performed without this constrain!');
	ErroR   = 1 ;
    end
end
fitness = originalFitness;
ParetoFrontier = zeros(ParetoFront(end),1);
for a = 1 : length(ParetoFrontier)   %% equal to number of Pareto front lines
    for b = 1 : length(ParetoFront)  %% equal to number of structures
	if ParetoFront(b) == a
	    ParetoFrontier(a) = ParetoFrontier(a) + 1;
	end
    end
end

WriteGenerationPareto(Ranking, ParetoFront, fitness, Enthalpies, Properties, paretoRanking, ORG_STRUC.opt_sign, ErroR);
%WriteParetoRanking(POP_STRUC.resFolder);

POP_STRUC.fitness = fitness;
POP_STRUC.ranking = Ranking;
POP_STRUC.paretoFront = ParetoFrontier;
POP_STRUC.bad_rank = bad_rank;
safesave ('Current_POP.mat', POP_STRUC)

%%%----------------------------------------------
function [Ranking, ParetoFront] = RANKING(BestRank, WorseRank)

if length(BestRank) ~= 0
   for a = size(WorseRank,1) : -1 : 1
     for b = 1 : size(BestRank,1)
        if WorseRank(a,1) == BestRank(b,1)
            WorseRank(a,:) = [];
            break;
        end
     end
   end
   if length(WorseRank) ~= 0
      WorseRank(:,2) = WorseRank(:,2) + BestRank(end,2);
      if WorseRank(1,2) ~= (BestRank(end,2) + 1)
         Diff = WorseRank(1,2) - (BestRank(end,2) + 1);
         WorseRank(:,2) = WorseRank(:,2) - Diff;
      end
      for a = 2 : size(WorseRank,1)
         if WorseRank(a,2) > WorseRank(a-1,2) + 1
            Diff = WorseRank(a,2) - (WorseRank(a-1,2) + 1);
            WorseRank(a:end,2) = WorseRank(a:end,2) - Diff;
         end
      end
   end
   result = [BestRank ; WorseRank];
   Ranking     = result(:,1);
   ParetoFront = result(:,2);
else
   Ranking     = WorseRank(:,1);
   ParetoFront = WorseRank(:,2);
end
