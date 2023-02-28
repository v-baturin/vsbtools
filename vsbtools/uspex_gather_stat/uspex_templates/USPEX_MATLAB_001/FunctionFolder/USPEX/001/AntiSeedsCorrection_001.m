function fitness = AntiSeedsCorrection_001(fitness)

global ORG_STRUC
global POP_STRUC
global ANTISEEDS

% AntiSeeds Policy 
%
%   antiSeedsActivation=N >0  : Gaussians are added to all structures starting from generation N
%   antiSeedsActivation=N <0  : Gaussians are only added to the best structure of each generation, starting from generation N
%   antiSeedsActivation=N =0  : Gaussians are only added to the structures put in the AntiSeeds folder

%finding antiSeedsMax
fitness_wk = []; %fitnesses without keptbest structures from previous generation
for i = 1 : length(fitness)
    if ~strcmp(POP_STRUC.POPULATION(i).howCome, 'keptBest') && ~strcmp(POP_STRUC.POPULATION(i).howCome, 'convexHull')
        fitness_wk(end+1) = fitness(i);
    end
end
number_kb = length(fitness) - length(fitness_wk); %number of keptBest and convexHull structures
fit_wk_len = round(length(fitness)*ORG_STRUC.bestFrac) - number_kb;
fitness_wk_sort = sort(fitness_wk);
antiSeedsMax = ORG_STRUC.antiSeedsMax*((sum(fitness_wk_sort(1:fit_wk_len))/fit_wk_len)-min(fitness)); 

fit_len = round(length(fitness)*ORG_STRUC.bestFrac);
% normalised on fitness variance (sort of)
%fitness_sort = sort(fitness);
%fitness_limit = fitness_sort(fit_len); %structures with fitness<=fitness_limit are important

% in this code we find ranking of the strutucres in current generation.
% The code repeats one from FitnessRankingCorrection function, but it is
% the simplest way not to change the other parts of teh code
[nothing, ranking] = sort(fitness);
Step = length(ORG_STRUC.abinitioCode);
[ranking, bad_rank, fitness1] = degradeSimilar(ORG_STRUC.toleranceFing, ranking, fitness, Step);

% antiSeedsMax = ORG_STRUC.antiSeedsMax*((sum(fitness_sort(1:fit_len))/fit_len)-min(fitness_sort)); 

% limit the minimal possible antiSeedsMax
if (POP_STRUC.generation > ORG_STRUC.antiSeedsActivation + 1) && ~isempty(ANTISEEDS)    
    if antiSeedsMax < 0.5*ANTISEEDS(end).Max0
        antiSeedsMax = 0.5*ANTISEEDS(end).Max0;
    end
    if antiSeedsMax > 2*ANTISEEDS(end).Max0
        antiSeedsMax = 2*ANTISEEDS(end).Max0;
    end    
end


count = 0;
antiSeedsSigma = 0;
for fit1 = 1 : fit_len - 1
    f1 = POP_STRUC.POPULATION(fit1).FINGERPRINT;
    numIons1 = POP_STRUC.POPULATION(fit1).numIons;
    
    for fit2 = fit1+1 : fit_len
        f2 = POP_STRUC.POPULATION(fit2).FINGERPRINT;
        numIons2 = POP_STRUC.POPULATION(fit2).numIons;        
        
        if ~isempty(f1) && ~isempty(f2) && isequal(numIons1, numIons2)
            weight = CalcWeight_001(numIons1);
            cos_dist = cosineDistance(f1, f2, weight);
            count = count + 1;
            antiSeedsSigma = antiSeedsSigma + cos_dist;
        end
    end
end
if (fit_len < 2) || (antiSeedsSigma < 0.000000000000001)
    antiSeedsSigma = 0.001;
    count = 1;
end
antiSeedsSigma = ORG_STRUC.antiSeedsSigma*antiSeedsSigma/count;


% initial the Sigma && Max at 1st genertion
if POP_STRUC.generation == 1
    for antiS_loop = 1 : length(ANTISEEDS)
        ANTISEEDS(antiS_loop).Sigma = antiSeedsSigma;
        ANTISEEDS(antiS_loop).Max   = antiSeedsMax;
        ANTISEEDS(antiS_loop).Max0  = antiSeedsMax; % without coefficient depending on composition
    end
end

if ~isempty(ANTISEEDS)
    if ~isempty(ANTISEEDS(1).FINGERPRINT) % ANTISEEDS exist
        for antiS_loop = 1 : length(ANTISEEDS)
            f1 = ANTISEEDS(antiS_loop).FINGERPRINT;
            numIons1 = ANTISEEDS(antiS_loop).numIons;
            
            for fit_loop = 1 : length(POP_STRUC.POPULATION)
                f2 = POP_STRUC.POPULATION(fit_loop).FINGERPRINT;
                numIons2 = POP_STRUC.POPULATION(fit_loop).numIons;
                
                if ~isempty(f1) && ~isempty(f2) && isequal(numIons1, numIons2)
                    weight = CalcWeight_001(numIons1);
                    cos_dist = cosineDistance(f1, f2, weight);
                    fitness(fit_loop) = fitness(fit_loop) + ANTISEEDS(antiS_loop).Max*exp(-cos_dist^2/(2*ANTISEEDS(antiS_loop).Sigma^2));
                end
            end
        end
    end
end

if ORG_STRUC.antiSeedsActivation > 0
    % add new antiseeds for ALL structures in the population
    if POP_STRUC.generation >= abs(ORG_STRUC.antiSeedsActivation)
        for i = 1 : length(POP_STRUC.POPULATION)
            f1 = POP_STRUC.POPULATION(i).FINGERPRINT;
            ANTISEEDS(end+1).FINGERPRINT = f1;
            ANTISEEDS(end).Sigma = antiSeedsSigma;
            ANTISEEDS(end).numIons = POP_STRUC.POPULATION(i).numIons;
            coef = find_coef_antiseeds(ANTISEEDS(end).numIons,fitness,ranking);
            ANTISEEDS(end).Max = coef * antiSeedsMax;
            ANTISEEDS(end).Max0 = antiSeedsMax;
        end
    end
elseif ORG_STRUC.antiSeedsActivation < 0
     % add new antiseed for the best structure only
     if POP_STRUC.generation >= abs(ORG_STRUC.antiSeedsActivation)
         [tmp,ind] = min(fitness);
         f1 = POP_STRUC.POPULATION(ind).FINGERPRINT;
         ANTISEEDS(end+1).FINGERPRINT = f1;
         ANTISEEDS(end).Sigma = antiSeedsSigma;
         ANTISEEDS(end).numIons = POP_STRUC.POPULATION(tmp).numIons;
         ANTISEEDS(end).Max = antiSeedsMax;
         ANTISEEDS(end).Max0 = antiSeedsMax;
         %status = ['AntiSeed added at generation' num2str(POP_STRUC.generation)]
     end
end

if ~isempty(ANTISEEDS)
    if isempty(ANTISEEDS(1).FINGERPRINT) % remove initial empty antiseed, if needed
        ANTISEEDS(1) = [];
    end
end

safesave ('ANTISEEDS.mat', ANTISEEDS)
