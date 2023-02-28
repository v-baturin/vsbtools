function createORG_MINHOP(inputFile)

global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

useBasicCell = python_uspex(getPy, ['-f ' inputFile ' -b retryBasicEvery -c 1']);
if ~isempty(useBasicCell)
   ORG_STRUC.useBasicCell = str2num(useBasicCell); % use basic cell for softmutations on top of the best structure, every 'ORG_STRUC.useBasicCell' generations
end

numGenerations = python_uspex(getPy, ['-f ' inputFile ' -b numGenerations -c 1']);
if isempty(numGenerations)
  numGenerations = '100'; % default
end
ORG_STRUC.numGenerations = str2num(numGenerations);

populationSize = python_uspex(getPy, ['-f ' inputFile ' -b populationSize -c 1']);
if isempty(populationSize)
  pS = 10*round(2*sum(ORG_STRUC.numIons)/10); % default
 if pS < 10
  pS = 10;
 elseif pS > 60
  pS = 60;
 end
 populationSize = num2str(pS);
end
%ORG_STRUC.populationSize = str2num(populationSize);
ORG_STRUC.populationSize = 2;

maxIncrease = python_uspex(getPy, ['-f ' inputFile ' -b maxCellIncrease -c 1']);
if isempty(maxIncrease)
% maxIncrease = ceil(ORG_STRUC.populationSize)^(1/3); % default, so that we can choose different supercell for every generation member
   ORG_STRUC.maxIncrease = 5; % default
else
   ORG_STRUC.maxIncrease = str2num(maxIncrease); % maximum cell multiplication factor when supercell is built
end

maxVectorLength = python_uspex(getPy, ['-f ' inputFile ' -b maxVectorLength -c 1']);
if isempty(maxVectorLength)
  lat_tmp = latConverter(ORG_STRUC.lattice);
  ORG_STRUC.maxVectorLength = 3*max(lat_tmp(1:3)); % default
else
  ORG_STRUC.maxVectorLength = str2num(maxVectorLength);
end

howManyMut = python_uspex(getPy, ['-f ' inputFile ' -b mutationDegree -c 1']);
if isempty(howManyMut)
  howManyMut = 0;
  for i = 1 : length(ORG_STRUC.atomType)
    howManyMut = howManyMut + str2num(covalentRadius(ORG_STRUC.atomType(i)));
  end
  ORG_STRUC.howManyMut = 3*howManyMut/length(ORG_STRUC.atomType);
else
  ORG_STRUC.howManyMut = str2num(howManyMut);
end

if ORG_STRUC.varcomp
  maxAt = python_uspex(getPy, ['-f ' inputFile ' -b maxAt -c 1']);
  if isempty(maxAt)
   maxAt = '100'; % default
  end
  ORG_STRUC.maxAt = str2num(maxAt);
else
  ORG_STRUC.maxAt = sum(ORG_STRUC.numIons);
end

minVectorLength = python_uspex(getPy, ['-f ' inputFile ' -b minVectorLength -c 1']);
if isempty(minVectorLength)
  minVectorLength = 0;
  for i = 1 : length(ORG_STRUC.atomType)
    if minVectorLength < 2*str2num(covalentRadius(ORG_STRUC.atomType(i)))
       minVectorLength = 2*str2num(covalentRadius(ORG_STRUC.atomType(i)));
    end
  end
  ORG_STRUC.minVectorLength = minVectorLength; % default
else
  ORG_STRUC.minVectorLength = str2num(minVectorLength);
end

useBasicCell = python_uspex(getPy, ['-f ' inputFile ' -b retryBasicEvery -c 1']);
if ~isempty(useBasicCell)
   ORG_STRUC.useBasicCell = str2num(useBasicCell); % use basic cell for softmutations on top of the best structure, every 'ORG_STRUC.useBasicCell' generations
end

