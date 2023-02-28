function update_STUFF(inputFile)

global ORG_STRUC
global POP_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

% READ FROM AUTOADAPT AND USERINPUT
numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if ~isempty(numParallelCalcs)
    ORG_STRUC.numParallelCalcs = str2num(numParallelCalcs);
end

% should the old guys be reoptimized?
reoptOld = python_uspex(getPy, ['-f ' inputFile ' -b reoptOld -c 1']);
if ~isempty(reoptOld)
    ORG_STRUC.reoptOld = str2num(reoptOld);
end

% maximum of gaussian for antiseed correction
antiSeedsMax = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsMax -c 1']);
if isempty(antiSeedsMax)
    antiSeedsMax = '0.000'; % default
end
ORG_STRUC.antiSeedsMax = str2num(antiSeedsMax);

% sigma of gaussian for antiseed correction
antiSeedsSigma = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsSigma -c 1']);
if isempty(antiSeedsSigma)
    antiSeedsSigma = '0.001'; % default
end
ORG_STRUC.antiSeedsSigma = str2num(antiSeedsSigma);

% switch the use of order in variation operators (except coormutation) on and off
ordering = python_uspex(getPy, ['-f ' inputFile ' -b ordering_active -c 1']);
if isempty(ordering)
    ordering = '1'; % default
end
ORG_STRUC.ordering = str2num(ordering);

% how many generations wait to add extra antiseed if the best structure is not changing
antiSeedsActivation = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsActivation -c 1']);
if isempty(antiSeedsActivation)
    antiSeedsActivation = '1'; % default
end
ORG_STRUC.antiSeedsActivation = str2num(antiSeedsActivation);

fracRand = python_uspex(getPy, ['-f ' inputFile ' -b fracRand -c 1']);
if isempty(fracRand)
    fracRand = '0.1'; % default
end
ORG_STRUC.fracRand = str2num(fracRand);

%%%%%%%%%%%%%%%%%%%% END defining tournement routine %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save ('Current_ORG.mat', 'ORG_STRUC')

