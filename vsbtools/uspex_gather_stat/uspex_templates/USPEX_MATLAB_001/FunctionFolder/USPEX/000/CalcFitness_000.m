function fitness = CalcFitness_000()

global POP_STRUC
global ORG_STRUC

fitness = zeros(1,length(POP_STRUC.POPULATION));
%Enthalpy
for fit_loop = 1:length(POP_STRUC.POPULATION)
    fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end);
end

fitness = ORG_STRUC.opt_sign*fitness; % change the mode of optimization (minimization <=> maximization)

fitnessRange = FitnessLimits(ORG_STRUC.opt_sign * ORG_STRUC.optType);
for i = 1 : length(fitness)
    if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999    % structure with errors
        fitness(i) = 100000;
    elseif fitness(i) < fitnessRange(1)
        fitness(i) = 100000;
    elseif fitness(i) > fitnessRange(2)
        fitness(i) = 100000;
    end
end

