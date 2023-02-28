function update_STUFF_old()
% $Rev$
% $Author$
% $Date$

global ORG_STRUC

popSize = ORG_STRUC.populationSize;

% sense. So in that case we stick to pure heredity.
if size(ORG_STRUC.numIons,2)==1
    ORG_STRUC.fracPerm = 0;
    ORG_STRUC.fracGene = ORG_STRUC.fracGene+ORG_STRUC.fracPerm;
end

howManyOffsprings     = round(popSize*ORG_STRUC.fracGene);
howManyPermutations   = round(popSize*ORG_STRUC.fracPerm);
howManyAtomMutations  = round(popSize*ORG_STRUC.fracAtomsMut);
howManyRand           = round(popSize*ORG_STRUC.fracRand);
howManyRotations      = round(popSize*ORG_STRUC.fracRotMut);      % molecules && proteins
howManyTransmutations = round(popSize*ORG_STRUC.fracTrans);       % varcom only
howManySpinmutations  = round(popSize*ORG_STRUC.fracSpin);       % spin operation
howManySecSwitch      = round(popSize*ORG_STRUC.fracSecSwitch);   % proteins only
howManyShiftBorder    = round(popSize*ORG_STRUC.fracShiftBorder); % proteins only

sum_offsprings = howManyOffsprings     + howManyPermutations  + ...
                 howManyTransmutations + howManyAtomMutations + ...
                 howManyRand           + howManyRotations     + ...
                 howManySpinmutations  + ...
                 howManySecSwitch      + howManyShiftBorder;

if sum_offsprings > popSize
    howManyOffsprings     = round((howManyOffsprings    *popSize)/sum_offsprings);
    howManyPermutations   = round((howManyPermutations  *popSize)/sum_offsprings);
    howManyAtomMutations  = round((howManyAtomMutations *popSize)/sum_offsprings);
    howManyRand           = round((howManyRand          *popSize)/sum_offsprings);
    howManyTransmutations = round((howManyTransmutations*popSize)/sum_offsprings);
    howManyRotations      = round((howManyRotations     *popSize)/sum_offsprings);
    howManySpinmutations  = round((howManySpinmutations *popSize)/sum_offsprings);
    howManySecSwitch      = round((howManySecSwitch     *popSize)/sum_offsprings);
    howManyShiftBorder    = round((howManyShiftBorder   *popSize)/sum_offsprings);
    sum_offsprings        = howManyOffsprings     + howManyPermutations  + ...
                            howManyTransmutations + howManyAtomMutations + ...
                            howManyRand           + howManyRotations     + ...
                            howManySpinmutations  + ...
                            howManySecSwitch      + howManyShiftBorder;
end

howManyleft = max([popSize-sum_offsprings, 0]);
if ~ORG_STRUC.constLattice
    howManyMutations = howManyleft;  %lat_Mutation
else
    howManyOffsprings = howManyOffsprings + howManyleft;
    howManyMutations = 0;  %lat_Mutation
end

ORG_STRUC.howManyMutations     = howManyMutations;
ORG_STRUC.howManyPermutations  = howManyPermutations;
ORG_STRUC.howManyAtomMutations = howManyAtomMutations;
ORG_STRUC.howManyRand          = howManyRand;
ORG_STRUC.howManyRotations     = howManyRotations;      % molecules && proteins
ORG_STRUC.howManyOffsprings    = howManyOffsprings;
ORG_STRUC.howManyTrans         = howManyTransmutations;
ORG_STRUC.howManySpinmutations = howManySpinmutations;% spin
ORG_STRUC.howManySecSwitch     = howManySecSwitch;      % proteins
ORG_STRUC.howManyShiftBorder   = howManyShiftBorder;    % proteins
ORG_STRUC.howManyRandTop       = 0;
ORG_STRUC.fracRandTop          = 0;
