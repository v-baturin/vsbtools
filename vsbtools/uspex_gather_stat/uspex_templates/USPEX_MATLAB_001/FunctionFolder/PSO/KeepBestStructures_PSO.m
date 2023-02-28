function addon_diff = KeepBestStructures_PSO()
%If you want to change how many steps you want to do for reoptimization,
%change this line
%OFF_STRUC.POPULATION(end).Step = length(ORG_STRUC.abinitioCode)-1;
% -1 means only relax for the last step

% this function lets best structures to survive, taking into acount dynamicalBestHM and other stuff
% 1-ORG_STRUC.bestFrac determines how many structures are thrown away

global ORG_STRUC
global POP_STRUC
global OFF_STRUC

keepBestHM = ORG_STRUC.keepBestHM;

if keepBestHM > round(ORG_STRUC.bestFrac*length(POP_STRUC.ranking))
    keepBestHM = round(ORG_STRUC.bestFrac*length(POP_STRUC.ranking));
end

% Add good guys from the old population and put their Step counter on the last Step
% We want only different guys to be added, if we use fingerprints

fitness = POP_STRUC.fitness;
first_mean = mean(fitness(POP_STRUC.ranking(1:round(ORG_STRUC.bestFrac*length(POP_STRUC.ranking)))));

% if we purely softmutate a generation we choose all different structures to survive
% if ORG_STRUC.dynamicalBestHM == 2 we have to tune tolerance so that exactly keepBestHM structures are chosen after clusterisation
% Tuning is done via binary search algorithm
tolerance = ORG_STRUC.toleranceBestHM;
deltaTol = tolerance/2;
decentRank = round(length(POP_STRUC.ranking)*ORG_STRUC.bestFrac);
doneOr = 0;
if ORG_STRUC.dynamicalBestHM == 2
    tolerance = 0.5;
    deltaTol = tolerance/2;
end

while doneOr == 0
    chosen = zeros(1,length(POP_STRUC.ranking));
    addon = 1;
    addon_diff = 0;
    while addon_diff < keepBestHM
        if addon > length(POP_STRUC.ranking)
            break;
        end
        if (addon > decentRank) & (ORG_STRUC.dynamicalBestHM > 0)
            break;
        end
        
        good_structure = 1;
        if (POP_STRUC.generation > ORG_STRUC.doFing-1)
            for j = 1 : length(POP_STRUC.ranking)
                if chosen(j) == 1
                    if (cosineDistance(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).FINGERPRINT, POP_STRUC.POPULATION(POP_STRUC.ranking(j)).FINGERPRINT, ORG_STRUC.weight) < tolerance)
                        good_structure = 0;
                    end
                end
            end
            if ((length(POP_STRUC.POPULATION)-addon)-(keepBestHM-addon_diff) <= 1) & (ORG_STRUC.dynamicalBestHM == 0)
                good_structure = 1;
            end
        end
        
        % dynamically vary keepBestHM
        if ORG_STRUC.dynamicalBestHM
            if addon > 1
                if fitness(POP_STRUC.ranking(addon)) > first_mean + power(var(fitness(POP_STRUC.ranking(1:round(ORG_STRUC.bestFrac*length(POP_STRUC.ranking))))), 0.5)
                    break;
                end;
            end;
        end
        
        if good_structure
            chosen(addon) = 1;
            addon_diff = addon_diff + 1;
        end
        addon = addon+1;
    end
    
    
    if ORG_STRUC.dynamicalBestHM < 2
        doneOr = 1;
    else
        if sum(chosen(:)) == keepBestHM
            doneOr = 1;
        elseif sum(chosen(:)) < keepBestHM
            tolerance = tolerance - deltaTol;
        else
            tolerance = tolerance + deltaTol;
        end
        deltaTol = deltaTol/2;
        if deltaTol < 0.000001
            doneOr = 1;
        end
    end
end

for addon = 1 : length(POP_STRUC.ranking)
    if chosen(addon)
        dummy = OFF_STRUC.POPULATION(end); % needed for octave compatibility
        OFF_STRUC.POPULATION(end+1) = dummy;
        OFF_STRUC.POPULATION(end)=POP_STRUC.POPULATION(POP_STRUC.ranking(addon));
        OFF_STRUC.POPULATION(end).Parents = [];
        OFF_STRUC.POPULATION(end).howCome = 'keptBest';
        
        info_parents = struct('parent',{}, 'enthalpy', {});
        info_parents(1).keptAsBest = num2str(POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).Number);
        info_parents.enthalpy = POP_STRUC.POPULATION(POP_STRUC.ranking(addon)).Enthalpies(end);
        OFF_STRUC.POPULATION(end).Parents = info_parents;
        if ORG_STRUC.reoptOld
            OFF_STRUC.POPULATION(end).Step = length(ORG_STRUC.abinitioCode)-1;
        else
            OFF_STRUC.POPULATION(end).Step = length(ORG_STRUC.abinitioCode)+1;
        end
    end
end

