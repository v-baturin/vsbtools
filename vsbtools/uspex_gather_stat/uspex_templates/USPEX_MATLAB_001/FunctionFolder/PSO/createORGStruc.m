function createORGStruc(inputFile)

% USPEX Version 9.4.1
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
warning off
global ORG_STRUC

ORG_STRUC = struct('dimension',{},'varcomp',{},'molecule',{},...                % EA
    'maxErrors',{},'numProcessors', {}, 'numParallelCalcs',{},...               % others
    'stopCrit',{}, 'reoptOld', {}, 'bestFrac',{},'repeatForStatistics',{},...
    'fracRand',{}, 'percSliceShift', {}, 'stopFitness', {},...
    'chargeNeutrality',{}, 'checkConnectivity',{}, 'howManyMut',{},...
    'valences',{}, 'NvalElectrons', {}, 'goodBonds',{}, ...
    'abinitioCode',{}, 'Kresol',{}, 'commandExecutable',{},...                  % abinitio calculation
    'remoteFolder', {}, 'platform',{}, 'ExternalPressure', {}, ...
    'latVolume',{},'constLat',{},'constLattice',{},...                          % lat and vol
    'splitInto', {}, ...                                                        % split cell
    'homePath', {}, 'USPEXPath',{}, 'specificFolder', {}, 'log_file', {}, ...   % path and output
    'atomType',{}, 'ionCh',{}, 'ionChange', {}, 'numIons', {},...               % atom
    'averageEnergy',{}, 'averageFitness',{},'correlation_coefficient',{},...    % aging
    'antiSeedsMax',{}, 'antiSeedsSigma',{}, 'ordering',{}, 'antiSeedsActivation',{},...
    'optType',{}, 'opt_sign',{},'cor_dir', {}, ...                              % optimization
    'RmaxFing',{}, 'deltaFing', {}, 'sigmaFing', {}, 'doFing', {},...           % Fingerprint
    'toleranceFing',{}, 'toleranceBestHM',{}, 'dynamicalBestHM',{}, ...
    'keepBestHM',{}, 'erf_table',{}, 'weight',{},...
    'PSO',{},'PSO_softMut',{},'PSO_BestStruc',{},'PSO_BestEver',{},...          % PSO
    'firstGeneSplit',{},'splitN',{},'minAt',{},'maxAt',{},'firstGeneMax',{},... % varcomp
    'doSpaceGroup',{}, 'nsym',{},'nsymN',{},...                                 % Symmetry
    'sym_coef',{}, 'symmetrize',{}, 'SGtolerance',{},...                        % Symmetry
    'pickUpYN',{},'pickUpGen',{}, 'pickUpFolder',{}, 'pickedUP', {},...                         % Pick up
    'minAngle',{}, 'minDiagAngle',{}, 'minDistMatrice',{},...                   % Constraints
    'minVectorLength',{},...                                                    % Constraints
    'numGenerations',{}, 'populationSize',{}, 'initialPopSize',{}, ...          % Pop and Generation
    'fixRndSeed',{});

% defaults
createORGDefault();

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%%%% Whether we fix the rand generator seeds
fixRndSeed = python_uspex(getPy, ['-f ' inputFile ' -b fixRandSeed -c 1']);
if ~isempty(fixRndSeed)
    ORG_STRUC.fixRndSeed = str2num(fixRndSeed);
end

%  [nothing, PSO_softMut] = unix(['./getStuff INPUT.txt PSO_softMut 1']);
PSO_softMut = python_uspex(getPy, ['-f ' inputFile ' -b PSO_softMut -c 1']);

% what % of structures are softmutated
ORG_STRUC.PSO_softMut = str2num(PSO_softMut);
%  [nothing, PSO_BestStruc] = unix (['./getStuff INPUT.txt PSO_BestStruc 1']);
PSO_BestStruc =python_uspex(getPy, ['-f ' inputFile ' -b PSO_BestStruc -c 1']);
% what % are breeded with best in their history
ORG_STRUC.PSO_BestStruc = str2num(PSO_BestStruc);
%  [nothing, PSO_BestEver] = unix (['./getStuff INPUT.txt PSO_BestEver 1']);
PSO_BestEver =python_uspex(getPy, ['-f ' inputFile ' -b PSO_BestEver -c 1']);
% what % are breeded with the very best structure
ORG_STRUC.PSO_BestEver = str2num(PSO_BestEver);
%%%%%%%%%%%%% EXPERIMENTAL - PSO %%%%%%%%%%%%%%%%%


% %%%%%%%% RESULTS
existfold = 1;
folderNum = 0;
while ~isempty(existfold)
    folderNum = folderNum + 1;
    ORG_STRUC.resFolder = ['results' num2str(folderNum)];
    [nothing1,existfold] = mkdir(ORG_STRUC.resFolder);
end
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'w');

% Create a header automatically:
description    = 'Particle Swarm Optimization Code for Structure Prediction';
funcFold       = [ORG_STRUC.USPEXPath '/FunctionFolder'];
formatted_rows = createHeader('USPEX_pic.txt', description, funcFold);

for i = 1 : size(formatted_rows, 2)
    fprintf(fp, [formatted_rows{i} '\n']);
end


% Cite:
text = {'Please cite the following suggested papers',       ...
    'when you publish the results obtained using', ...
    'Corrected PSO method within USPEX::' ...
    };
formatted_rows = createHeader_wrap(text);

text = {'Corrected PSO Method'};
formatted_rows = createHeader_wrap(text);
for i = 1 : size(formatted_rows, 2)
    fprintf(fp, [formatted_rows{i} '\n']);
end

text = {'Lyakhov A.O., Oganov A.R., Stokes H., Zhu Q., (2013)', ...
    'New developments in evolutionary structure prediction algorithm USPEX.',  ...
    'Comp. Phys. Comm., 184, 1172-1182' ...
    };
formatted_rows = createHeader_wrap(text, 'left');
for i = 1 : size(formatted_rows, 2)
    fprintf(fp, [formatted_rows{i} '\n']);
end

fprintf(fp, '\n');

fclose(fp);

%%%%%%% GET the code and cluster.
%[nothing, abinitioCode] = unix(['./getFromTo abinitioCode ENDabinit ' inputFile]);
%ORG_STRUC.abinitioCode = str2num(abinitioCode);
ORG_STRUC.abinitioCode = python_uspex(getPy, ['-f ' inputFile ' -b abinitioCode -e ENDabinit'], 1);
ExternalPressure = python_uspex(getPy, ['-f ' inputFile ' -b ExternalPressure -c 1']);
if ~isempty(ExternalPressure)
    ORG_STRUC.ExternalPressure = str2num(ExternalPressure);
end


%[nothing, whichCluster] = unix(['./getStuff ' inputFile ' whichCluster 1']);
whichCluster = python_uspex(getPy, ['-f ' inputFile ' -b whichCluster -c 1']);
if ~isempty(whichCluster)
    if ~isempty(str2num(whichCluster)) % number
        ORG_STRUC.platform = abs(str2num(whichCluster));
    else
        whichCluster(end) = [];
        if strcmp(whichCluster, 'nonParallel')
            ORG_STRUC.platform = 0;
        elseif strcmp(whichCluster, 'CFN')
            ORG_STRUC.platform = 3;
        elseif strcmp(whichCluster, 'QSH')
            ORG_STRUC.platform = 4;
        elseif strcmp(whichCluster, 'QSH2')
            ORG_STRUC.platform = 5;
        elseif strcmp(whichCluster, 'xservDE')
            ORG_STRUC.platform = 6;
        elseif strcmp(whichCluster, 'MIPT')
            ORG_STRUC.platform = 7;
        end
    end
end
if ORG_STRUC.platform > 7
    disp(['you did not chose a valid platform. Program STOPS...']);
    fprintf(fp,'you did not chose a valid platform. Program STOPS...\n');
    quit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% optimization parameter; optimize by: enthalpy, volume, hardness, etc
optType = python_uspex(getPy, ['-f ' inputFile ' -b optType -c 1']);
if ~isempty(optType)
    if optType(1) == '-'  % non documented: '-' at start means 'reverse' optimization (minimization <=> maximization)
        ORG_STRUC.opt_sign = -1;
    end
    optType(end) = [];
    
    if ~isempty(str2num(optType)) % number
        ORG_STRUC.optType=abs(str2num(optType));
    else
        if strcmpi(optType, 'enthalpy')
            ORG_STRUC.optType = 1;
        elseif strcmpi(optType, 'volume')
            ORG_STRUC.optType = 2;
        elseif strcmpi(optType, 'hardness')
            ORG_STRUC.optType = 3;
        elseif strcmpi(optType, 'struc_order')
            ORG_STRUC.optType = 4;
        elseif strcmpi(optType, 'aver_dist')
            ORG_STRUC.optType = 5;
        elseif strcmpi(optType, 'diel_sus')
            ORG_STRUC.optType = 6;
        elseif strcmpi(optType, 'gap')
            ORG_STRUC.optType = 7;
        elseif strcmpi(optType, 'diel_gap')
            ORG_STRUC.optType = 8;
        elseif strcmpi(optType, 'mag_moment')
            ORG_STRUC.optType = 9;
        elseif strcmpi(optType, 'quasientropy')
            ORG_STRUC.optType = 10;
        end
    end
    
    optType_OK = [1 : 10, 1101 : 1111];
    if isempty( find( ORG_STRUC.optType==optType_OK) )
        fprintf(fp, 'you did not chose a valid optType. Program STOPS...\n');
        disp(['you did not chose a valid optType. Program STOPS...']);
        quit;
    end
end

% 1:'enthalpy'
% 2:'volume'
% 3:'hardness'
% 4:'struc_order'
% 5:'aver_dist'
% 6:'diel_sus'
% 7:'gap'
% 8:'diel_gap'
% 9:'mag_moment'
%10:'quasientropy'

%%%%%CalcType
%3 digits (dimension: 0-3; molecule: 0/1; varcomp: 0/1)
calculationType = python_uspex(getPy, ['-f ' inputFile ' -b calculationType -c 1']);
if ~isempty(calculationType)
    if calculationType(1) == '-'  % 2D
        calculationType(1) = [];
        ORG_STRUC.dimension = -str2num(calculationType(1));
    else
        ORG_STRUC.dimension = str2num(calculationType(1));
    end
    
    ORG_STRUC.molecule  = str2num(calculationType(2));
    ORG_STRUC.varcomp   = str2num(calculationType(3));
end
%[nothing, howManySliceIterations] = unix (['./getStuff ' inputFile ' howManySliceIterations 1']);
howManySliceIterations = python_uspex(getPy, ['-f ' inputFile ' -b howManySliceIterations -c 1']);
if ~isempty(howManySliceIterations)
    ORG_STRUC.Nslabs = str2num(howManySliceIterations);
end
%[nothing, maxDistHeredity] = unix (['./getStuff ' inputFile ' maxDistHeredity 1']);
maxDistHeredity = python_uspex(getPy, ['-f ' inputFile ' -b maxDistHeredity -c 1']);
if ~isempty(maxDistHeredity)
    ORG_STRUC.maxDistHeredity = str2num(maxDistHeredity);
end
%[nothing, doSpaceGroup] = unix (['./getStuff ' inputFile ' doSpaceGroup 1']);
doSpaceGroup = python_uspex(getPy, ['-f ' inputFile ' -b doSpaceGroup -c 1']);
if isempty(doSpaceGroup)
    if (ORG_STRUC.dimension==0) | (ORG_STRUC.dimension==2) | (ORG_STRUC.dimension==-3)
        doSpaceGroup = '0'; % default
    else
        doSpaceGroup = '1'; % default
    end
end
ORG_STRUC.doSpaceGroup = str2num(doSpaceGroup);

% number of atoms to split cell to
%[nothing, splitInto] = unix(['./getFromTo splitInto EndSplitInto ' inputFile]);
splitInto =  python_uspex(getPy, ['-f ' inputFile ' -b splitInto -e EndSplitInto']);
if ~isempty(splitInto)
    ORG_STRUC.splitInto = str2num(splitInto);
end

% how many times to repeat the calculations with the same initial conditions
% this option allows to make some statistics rerunning USPEX a few times after calculations are done
repeatForStatistics = python_uspex(getPy, ['-f ' inputFile ' -b repeatForStatistics -c 1']);
if ~isempty(repeatForStatistics)
    ORG_STRUC.repeatForStatistics = str2num(repeatForStatistics);
end
repeatForStatistics = num2str(ORG_STRUC.repeatForStatistics);
if ORG_STRUC.repeatForStatistics > 1
    if ~exist('multiple_runs')
        unixCmd(['echo "' repeatForStatistics ' : runs_left" > multiple_runs']);
    end
end

%% fingerprints:
%[nothing, RmaxFing] = unix (['./getStuff ' inputFile ' RmaxFing 1']);
RmaxFing =  python_uspex(getPy, ['-f ' inputFile ' -b RmaxFing -c 1']);
if ~isempty(RmaxFing)
    ORG_STRUC.RmaxFing = str2num(RmaxFing);
end
if ORG_STRUC.RmaxFing > 0
    ORG_STRUC.doFing = 1;
else
    ORG_STRUC.doFing = 0;
end

%[nothing, deltaFing] = unix (['./getStuff ' inputFile ' deltaFing 1']);
deltaFing = python_uspex(getPy, ['-f ' inputFile ' -b deltaFing -c 1']);
if ~isempty(deltaFing)
    ORG_STRUC.deltaFing = str2num(deltaFing);
end

%[nothing, sigmaFing] = unix (['./getStuff ' inputFile ' sigmaFing 1']);
sigmaFing = python_uspex(getPy, ['-f ' inputFile ' -b sigmaFing -c 1']);
if ~isempty(sigmaFing)
    ORG_STRUC.sigmaFing = str2num(sigmaFing);
end

%[nothing, toleranceFing] = unix (['./getStuff ' inputFile ' toleranceFing 1']);
toleranceFing = python_uspex(getPy, ['-f ' inputFile ' -b toleranceFing -c 1']);
if ~isempty(toleranceFing)
    ORG_STRUC.toleranceFing = str2num(toleranceFing);
end

toleranceBestHM = python_uspex(getPy, ['-f ' inputFile ' -b toleranceBestHM -c 1']);
if ~isempty(toleranceBestHM)
    ORG_STRUC.toleranceBestHM = str2num(toleranceBestHM);
end
% maximum of gaussian for antiseed correction
%[nothing, antiSeedsMax] = unix (['./getStuff ' inputFile ' antiSeedsMax 1']);
antiSeedsMax = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsMax -c 1']);
if ~isempty(antiSeedsMax)
    ORG_STRUC.antiSeedsMax = str2num(antiSeedsMax);
end
% sigma of gaussian for antiseed correction
%[nothing, antiSeedsSigma] = unix (['./getStuff ' inputFile ' antiSeedsSigma 1']);
antiSeedsSigma = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsSigma -c 1']);
if ~isempty(antiSeedsSigma)
    ORG_STRUC.antiSeedsSigma = str2num(antiSeedsSigma);
end
% switch the use of order in variation operators (except coormutation) on and off
%[nothing, ordering] = unix (['./getStuff ' inputFile ' ordering_active 1']);
ordering = python_uspex(getPy, ['-f ' inputFile ' -b ordering_active -c 1']);
if ~isempty(ordering)
    ORG_STRUC.ordering = str2num(ordering);
end
% how many generations wait to add extra antiseed if the best structure is not changing
%[nothing, antiSeedsActivation] = unix (['./getStuff ' inputFile ' antiSeedsActivation 1']);
antiSeedsActivation = python_uspex(getPy, ['-f ' inputFile ' -b antiSeedsActivation -c 1']);
if ~isempty(antiSeedsActivation)
    ORG_STRUC.antiSeedsActivation = str2num(antiSeedsActivation);
end

% do we want to symmetrize the structure using Stokes SG determination code (using symmetrized.cif)
% done at the last optimization step (when not shaking - symmetrize)
symmetrize = python_uspex(getPy, ['-f ' inputFile ' -b symmetrize -c 1']);
if ~isempty(symmetrize)
    ORG_STRUC.symmetrize = str2num(symmetrize);
end
% Space group determination tolerance
%[noathing, SGtolerance] = unix(['./getStuff ' inputFile ' SymTolerance 1']);
SGtolerance = python_uspex(getPy, ['-f ' inputFile ' -b SymTolerance -c 1']);
if ~isempty(SGtolerance)
    SGtolerance = deblank(SGtolerance);
    if strcmp(lower(SGtolerance), 'high')
        ORG_STRUC.SGtolerance = 0.05;
    elseif strcmp(lower(SGtolerance), 'medium')
        ORG_STRUC.SGtolerance = 0.10;
    elseif strcmp(lower(SGtolerance), 'low')
        ORG_STRUC.SGtolerance = 0.20;
    else % number!
        ORG_STRUC.SGtolerance = str2num(SGtolerance);
    end
end
if isempty(ORG_STRUC.SGtolerance) %In case some illegal strings
    ORG_STRUC.SGtolerance = 0.10;
end

% what point symmetries have to be satisfied by clusters created randomly
% for crystals - what symmetry groups the randomly created crystals should belong to, format : 1-3 5 6-9 15  etc
nsym = python_uspex(getPy, ['-f ' inputFile ' -b symmetries -e endSymmetries']);
if isempty(nsym)
    if ORG_STRUC.dimension == 0
        nsym = 'E C2 D2 C4 C3 C6 T S2 Ch1 Cv2 S4 S6 Ch3 Th Ch2 Dh2 Ch4 D3 Ch6 O D4 Cv3 D6 Td Cv4 Dd3 Cv6 Oh Dd2 Dh3 Dh4 Dh6 Oh C5 S5 S10 Cv5 Ch5 D5 Dd5 Dh5 I Ih '; % default
    elseif ORG_STRUC.dimension == 3
        nsym = '2-230'; % default
    elseif ORG_STRUC.dimension == -2 | ORG_STRUC.dimension == 1
        nsym = '2-17';
    end
end
if ORG_STRUC.dimension == 0
    nsym(end) = [];
    ORG_STRUC.nsym = nsym;
    c1 = findstr(nsym, ' ');
    c = sort(str2num(['0 ' num2str(c1)]));
    c(end+1) = length(nsym) + 1;
    ind = zeros(1,2);
    indN = 0;
    for i = 2 : length(c)
        if c(i-1) + 1 > c(i) - 1
            continue
        end
        indN = indN + 1;
        ind(1, 1) = c(i-1) + 1;
        ind(1, 2) = c(i) - 1;
        if indN == 1
            ORG_STRUC.nsymN = ind;
        else
            ORG_STRUC.nsymN = cat(1, ORG_STRUC.nsymN, ind);
        end
    end
    if isempty(ORG_STRUC.nsym)
        ORG_STRUC.nsym = 'E';
        ORG_STRUC.nsymN = [1 1];
    end
elseif ORG_STRUC.dimension ~= 2
    ORG_STRUC.nsym = zeros(1, 230);
    c1 = findstr(nsym, ' ');
    c2 = findstr(nsym, '-');
    c = sort(str2num(['0 ' num2str(c1) ' ' num2str(c2)]));
    c(end+1) = length(nsym) + 1;
    ind1 = 1;
    for i = 2 : length(c)
        if c(i-1)+1 > c(i)-1
            continue
        end
        ind2 = str2num(nsym(c(i-1)+1 : c(i)-1));
        if ind2 == 0
            ind1 = 1;
            continue
        end
        if ~isempty(find(c2 == c(i-1)))
            for j = ind1 : ind2
                ORG_STRUC.nsym(j) = 1;
            end
        else
            ORG_STRUC.nsym(ind2) = 1;
        end
        ind1 = ind2;
    end
    if sum(ORG_STRUC.nsym) == 0
        ORG_STRUC.nsym(1) = 1;
    end
    ORG_STRUC.nsymN = [0 0];
end
% coefficient between mindist and symmetrization distance (1 by default, sometimes > 1 needed)
sym_coef = python_uspex(getPy, ['-f ' inputFile ' -b constraint_enhancement -c 1']);
if ~isempty(sym_coef)
    ORG_STRUC.sym_coef = str2num(sym_coef);
end
% number of ions of each type
%[nothing, numIons] = unix(['./getFromTo numSpeci EndNumSpeci ' inputFile]);
%ORG_STRUC.numIons = str2num(numIons);
numIons = python_uspex(getPy, ['-f ' inputFile ' -b numSpeci -e EndNumSpeci'], 1);
ORG_STRUC.numIons = numIons;

% types of ions; can be number, short name of full name. USPEX will use the numbers only
%[nothing, atomType] = unix(['./getFromTo atomType EndAtomType ' inputFile]);
atomType = python_uspex(getPy, ['-f ' inputFile ' -b atomType -e EndAtomType']);
atomType1 = zeros(1,size(ORG_STRUC.numIons,2));
c1 = findstr(atomType, ' ');
c = sort(str2num(['0 ' num2str(c1)]));
c(end+1) = length(atomType) + 1;
ind1 = 1;
for i = 2 : length(c)
    if c(i-1)+1 > c(i)-1
        continue
    end
    tmp = atomType(c(i-1)+1 : c(i)-1);
    tmp = strtrim(tmp);
    if ~isempty(str2num(tmp)) % number
        atomType1(ind1) = str2num(tmp);
    else
        if strcmp(tmp,'H.5')
            atomType1(ind1)=0.5;
        elseif strcmp(tmp,'H.75')
            atomType1(ind1)=0.75;
        elseif strcmp(tmp,'H1.25')
            atomType1(ind1)=1.25;
        elseif strcmp(tmp,'H1.5')
            atomType1(ind1)=1.5;
        else
            for j = 1 : 105
                if strcmp(lower(tmp), lower(elementFullName(j))) | strcmp(lower(tmp), lower(megaDoof(j)))
                    atomType1(ind1) = j;
                    break;
                end
            end
        end
    end
    ind1 = ind1 + 1;
end

for i=1:length(atomType1)
    if atomType1(i)==0
        atomType1(i)=[];
    end
end
ORG_STRUC.atomType = atomType1;

% valences for each type of atoms
valences = python_uspex(getPy, ['-f ' inputFile ' -b valences -e endValences']);
ORG_STRUC.valences = str2num(valences);
if isempty(valences) % default
    ORG_STRUC.valences = zeros(1,length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        if ORG_STRUC.atomType(i) == 0.5
            ORG_STRUC.valences(i) = 0.5;
        elseif ORG_STRUC.atomType(i) == 0.75
            ORG_STRUC.valences(i) = 0.75
        elseif ORG_STRUC.atomType(i) == 1.25
            ORG_STRUC.valences(i) = 1.25;
        elseif ORG_STRUC.atomType(i) == 1.5
            ORG_STRUC.valences(i) = 1.5;
        else
            ORG_STRUC.valences(i) = str2num(valence(ORG_STRUC.atomType(i)));
        end
    end
else
    ORG_STRUC.valences = str2num(valences);
end

% number of valence electrons for each type of atoms
NvalElectrons = python_uspex(getPy, ['-f ' inputFile ' -b valenceElectr -e endValenceElectr']);
if isempty(NvalElectrons) % default
    ORG_STRUC.NvalElectrons = zeros(1,length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        if ORG_STRUC.atomType(i) == 0.5
            ORG_STRUC.NvalElectrons(i) = 0.5;
        elseif ORG_STRUC.atomType(i) == 0.75
            ORG_STRUC.NvalElectrons(i) = 0.75
        elseif ORG_STRUC.atomType(i) == 1.25
            ORG_STRUC.NvalElectrons(i) = 1.25;
        elseif ORG_STRUC.atomType(i) == 1.5
            ORG_STRUC.NvalElectrons(i) = 1.5;
        else
            ORG_STRUC.NvalElectrons(i) = str2num(valenceElectronsNumber(ORG_STRUC.atomType(i)));
        end
    end
else
    ORG_STRUC.NvalElectrons = str2num(NvalElectrons);
end

% how many individuals per generation
populationSize = python_uspex(getPy, ['-f ' inputFile ' -b populationSize -c 1']);
if isempty(populationSize)
    if ORG_STRUC.varcomp == 1
        pS = 10 * round(2 * ORG_STRUC.maxAt / 10); % default
    else
        pS = 10 * round(2 * sum(ORG_STRUC.numIons) / 10); % default
    end
    if pS < 10
        pS = 10;
    elseif pS > 60
        pS = 60;
    end
    populationSize = num2str(pS);
end
ORG_STRUC.populationSize = str2num(populationSize);

% how many individuals in the first generation
initialPopSize = python_uspex(getPy, ['-f ' inputFile ' -b initialPopSize -c 1']);
if isempty(initialPopSize)
    initialPopSize = num2str(ORG_STRUC.populationSize); % default
end
ORG_STRUC.initialPopSize = str2num(initialPopSize);

maxErrors = python_uspex(getPy, ['-f ' inputFile ' -b maxErrors -c 1']);
if ~isempty(maxErrors)
    ORG_STRUC.maxErrors = str2num(maxErrors);
end
% connectivity check in hardness, distanceCheck etc
checkConnectivity = python_uspex(getPy, ['-f ' inputFile ' -b checkConnectivity -c 1']);
if isempty(checkConnectivity)
    if ORG_STRUC.optType == 3
        ORG_STRUC.checkConnectivity = 1;  %default
    else
        ORG_STRUC.checkConnectivity = 0;  %default
    end
else
    ORG_STRUC.checkConnectivity = str2num(checkConnectivity);
end
% how big a fraction of a generation shall be included
% in the process of heredity (the best ones of course)

fracRand = python_uspex(getPy, ['-f ' inputFile ' -b fracRand -c 1']);
if isempty(fracRand)
    ORG_STRUC.fracRand = str2num(fracRand);
end

howManyMut = python_uspex(getPy, ['-f ' inputFile ' -b mutationDegree -c 1']);
if isempty(howManyMut)
    howManyMut = 0;
    for i = 1 : length(ORG_STRUC.atomType)
        howManyMut = howManyMut + str2num(covalentRadius(ceil(ORG_STRUC.atomType(i))));
    end
    ORG_STRUC.howManyMut = 3 * howManyMut / length(ORG_STRUC.atomType);
else
    ORG_STRUC.howManyMut = str2num(howManyMut);
end

% how many generations shall be calculated
numGenerations = python_uspex(getPy, ['-f ' inputFile ' -b numGenerations -c 1']);
if ~isempty(numGenerations)
    ORG_STRUC.numGenerations = str2num(numGenerations);
end

% how big a fraction of structures generated by heredity shall be shifted in all dimensions
percSliceShift = python_uspex(getPy, ['-f ' inputFile ' -b percSliceShift -c 1']);
if ~isempty(percSliceShift)
    ORG_STRUC.percSliceShift = str2num(percSliceShift);
end
% number of repeating best structure that stops the calculation
stopCrit = python_uspex(getPy, ['-f ' inputFile ' -b stopCrit -c 1']);
if isempty(stopCrit)
    if ORG_STRUC.varcomp
        ORG_STRUC.stopCrit = ORG_STRUC.maxAt; % default
    else
        ORG_STRUC.stopCrit = sum(ORG_STRUC.numIons); % default
    end
else
    ORG_STRUC.stopCrit = str2num(stopCrit);
end

% information on the lattice
ORG_STRUC.lattice = python_uspex(getPy, ['-f ' inputFile ' -b Latticevalues -e Endvalues'], 1);

minVectorLength = python_uspex(getPy, ['-f ' inputFile ' -b minVectorLength -c 1']);
if isempty(minVectorLength)
    minVectorLength = 0;
    Vector = zeros(1, length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        Vector(i) = 2 * str2num(covalentRadius(ceil(ORG_STRUC.atomType(i))));
    end
    if ORG_STRUC.varcomp == 0
        minVectorLength = max(Vector);
    else
        minVectorLength = min(Vector);
    end
    ORG_STRUC.minVectorLength = minVectorLength; % default
else
    ORG_STRUC.minVectorLength = str2num(minVectorLength);
end

try
    specificSwaps = python_uspex(getPy, ['-f ' inputFile ' -b specificSwaps -e EndSpecific']);
    ORG_STRUC.specificSwaps = str2num(specificSwaps);
    if isempty(ORG_STRUC.specificSwaps)
        ORG_STRUC.specificSwaps = 0;
    end
catch
    ORG_STRUC.specificSwaps = 0;
end

% how many[nothing, keepBestHM] = unix (['./getStuff ' inputFile ' keepBestHM 1']); best guys to be taken in account for new pop
keepBestHM =  python_uspex(getPy, ['-f ' inputFile ' -b keepBestHM -c 1']);
if isempty(keepBestHM)
    keepBestHM = num2str(round(0.15 * ORG_STRUC.populationSize)); % default
end
ORG_STRUC.keepBestHM = str2num(keepBestHM);
if ORG_STRUC.keepBestHM < 1
    ORG_STRUC.keepBestHM = 1;
end

% dynamicalBestHM = 1 means that value of keepBestHM can vary during calculations with initial ORG_STRUC.keepBestHM as upper bound
% = 2 - clusterization algorithm used to keep exactly keepBestHM structures, as diverse as possible
dynamicalBestHM = python_uspex(getPy, ['-f ' inputFile ' -b dynamicalBestHM -c 1']);
if ~isempty(dynamicalBestHM)
    ORG_STRUC.dynamicalBestHM = str2num(dynamicalBestHM);
end
% should the old guys be reoptimized?
reoptOld = python_uspex(getPy, ['-f ' inputFile ' -b reoptOld -c 1']);
if ~isempty(reoptOld)
    ORG_STRUC.reoptOld = str2num(reoptOld);
end
% minimal possible distances between ions.
hardCore = python_uspex(getPy, ['-f ' inputFile ' -b IonDistances -e EndDistances']);
if isempty(hardCore)
    hardCore = zeros(1, length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        if ORG_STRUC.atomType(i) < 2.5
            hardCore = 0.3;
        else
            s = covalentRadius(ceil(ORG_STRUC.atomType(i)));
            hardCore(i) = str2num(s) * 0.6; %to set the possible ionDistance
        end
    end
else
    hardCore = str2num(hardCore);
end


%%%% check whether hard core radii or distance matrix were given and create matrix if necessary %%%%%%%
if size(hardCore,1) == size(hardCore,2)
    for i = 1 : length(ORG_STRUC.atomType)
        for j = i : length(ORG_STRUC.atomType)
            ORG_STRUC.minDistMatrice(i,j) = hardCore(i,j);
            ORG_STRUC.minDistMatrice(j,i) = hardCore(i,j);
        end
    end
else
    ORG_STRUC.minDistMatrice = zeros(length(hardCore));
    for hardInd_1 = 1 : length(hardCore)
        for hardInd_2 = hardInd_1 : length(hardCore)
            ORG_STRUC.minDistMatrice(hardInd_1,hardInd_2) = hardCore(hardInd_1) + hardCore(hardInd_2);
            ORG_STRUC.minDistMatrice(hardInd_2,hardInd_1) = hardCore(hardInd_1) + hardCore(hardInd_2);
        end
    end
end
%%%%%%%%%%%%%% END hard core %%%%%%%%%%%%%%%%%%

% bonds with nu higher than this value are always included into hardness formula
ORG_STRUC.goodBonds = zeros(length(ORG_STRUC.atomType));
goodBonds = python_uspex(getPy, ['-f ' inputFile ' -b goodBonds -e EndGoodBonds']);
if isempty(goodBonds)
    goodBonds = '0.15'; % default
end
gB = str2num(goodBonds);
if size(gB, 1) == 1        % just 1 value present or 1 type of atoms
    for i = 1 : length(ORG_STRUC.atomType)
        for j = 1 : length(ORG_STRUC.atomType)
            ORG_STRUC.goodBonds(i,j) = gB;
        end
    end
else
    for i = 1 : length(ORG_STRUC.atomType)
        for j = i : length(ORG_STRUC.atomType)
            ORG_STRUC.goodBonds(i,j) = gB(i,j);
            ORG_STRUC.goodBonds(j,i) = gB(i,j);
        end
    end
end

% Kresol is the desired resolution for KPOINTS
Kresol = python_uspex(getPy, ['-f ' inputFile ' -b KresolStart -e Kresolend']);
if isempty(Kresol) % default
    Ltmp = length(ORG_STRUC.abinitioCode);
    ORG_STRUC.Kresol = 0.2 * ones(1, Ltmp);
    if Ltmp > 1
        for i = 1 : Ltmp
            ORG_STRUC.Kresol(i) = 0.2 - (i-1) * 0.12 / (Ltmp-1);
        end
    end
else
    ORG_STRUC.Kresol = str2num(Kresol);
end


% if pickUpYN~=0 , then a previous calculation will be continued
pickUpYN = python_uspex(getPy, ['-f ' inputFile ' -b pickUpYN -c 1']);
if ~isempty(pickUpYN)
    ORG_STRUC.pickUpYN = str2num(pickUpYN);
end
% from which generation shall the last calculation be continued? If = 0, then the last available generation is taken
pickUpGen = python_uspex(getPy, ['-f ' inputFile ' -b pickUpGen -c 1']);
if ~isempty(pickUpGen)
    ORG_STRUC.pickUpGen = str2num(pickUpGen);
end
if ORG_STRUC.pickUpYN == 0
    ORG_STRUC.pickUpGen = 1; % needed for BestEnthalpies plot
end

% number of the results folder to be used. If = 0 , then the highest
% existing number is taken
pickUpFolder = python_uspex(getPy, ['-f ' inputFile ' -b pickUpFolder -c 1']);
if ~isempty(pickUpFolder)
    ORG_STRUC.pickUpFolder = str2num(pickUpFolder);
end

remoteFolder1 = python_uspex(getPy, ['-f ' inputFile ' -b remoteFolder -c 1']);
ORG_STRUC.remoteFolder = remoteFolder1(1 : length(remoteFolder1)-1);

% number of processors involved in each calculation
numProcessors = python_uspex(getPy, ['-f ' inputFile ' -b numProcessors -e EndProcessors']);
if isempty(numProcessors)
    ORG_STRUC.numProcessors = 8 * ones(1, length(ORG_STRUC.abinitioCode)); % default
else
    ORG_STRUC.numProcessors = str2num(numProcessors);
end

numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if ~isempty(numParallelCalcs)
    ORG_STRUC.numParallelCalcs = str2num(numParallelCalcs);
end

%%%%%%%%% END get all parameters from text file %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%% Command stuff
if ORG_STRUC.platform==0
    
    commandExecutable = python_uspex(getPy, ['-f ' inputFile ' -b commandExecutable -e EndExecutable']);
    
    helper = double(commandExecutable);
    tempIND = find (helper == 10);
    helpIND = zeros(length(tempIND)+1, 1);
    helpIND(2: end) = tempIND';
    for listLoop = 1:length(helpIND)-1
        listCommand{listLoop} = commandExecutable(helpIND(listLoop)+1 : helpIND(listLoop+1)-1);
    end
    
    for commandLoop = 1: length(ORG_STRUC.abinitioCode)
        if ORG_STRUC.abinitioCode(commandLoop) == 0
            tempCom{commandLoop} = '0';     % means we do not do any optimization at all! (used in order optimization)
        elseif listLoop < length(ORG_STRUC.abinitioCode)
            tempCom{commandLoop} = listCommand{1};
        else
            tempCom{commandLoop} = listCommand{commandLoop};
        end
    end
    ORG_STRUC.commandExecutable = tempCom;
end
%%%%%%%%%%%%% END command stuff

%%%%%%% LATTICE manipulations %%%%%%%%%%%%
% check the rare case for varcomp with 6 atomic types
varcomp_6_fixlat = 1;
if (ORG_STRUC.varcomp == 1) & (length(ORG_STRUC.atomType) == 6)
    if max(ORG_STRUC.lattice(4:6)) < 30   % 6 atomic volumes. If any value > 30 then it's fixed lattice instead
        varcomp_6_fixlat = 0;
    end
end

if size(ORG_STRUC.lattice, 1) == 1
    if (size(ORG_STRUC.lattice, 2) == 6) & (varcomp_6_fixlat == 1)    % fixed lattice in a form {a b c alpha beta gamma}?
        ORG_STRUC.lattice(4:6) = pi / 180 * ORG_STRUC.lattice(4 : 6);         % Matlab works with radians, input - in degrees
        ORG_STRUC.latVolume = det(latConverter(ORG_STRUC.lattice));
        ORG_STRUC.constLattice = 1;
        ORG_STRUC.constLat = latConverter(ORG_STRUC.lattice);   % we like 3x3 form of the lattice
        ORG_STRUC.lattice = []; % Octave compatibility
        ORG_STRUC.lattice = ORG_STRUC.constLat;   % 3x3 form
    else
        ORG_STRUC.latVolume = ORG_STRUC.lattice;
        ORG_STRUC.constLattice = 0;
    end
elseif sum(size(ORG_STRUC.lattice)) == 6
    ORG_STRUC.latVolume = det(ORG_STRUC.lattice);
    ORG_STRUC.constLattice = 1;
    ORG_STRUC.constLat = ORG_STRUC.lattice;
elseif ORG_STRUC.dimension ~= 2
    disp(['Please Specify the lattice values in INPUT.txt']);
end
%%%%%%% END LATTICE manipulations %%%%%%%%%%%%

if ORG_STRUC.varcomp == 1 | ORG_STRUC.dimension == 2 %varcomp and surface
    ORG_STRUC.weight = 1;
else
    if ORG_STRUC.molecule == 1
        numIons = ORG_STRUC.numMols;
    else
        numIons = ORG_STRUC.numIons;
    end
    % weight needed for normalisation of the cosine distance between fingerprints
    L = length(numIons);
    S = 0;
    ORG_STRUC.weight = zeros(L*L, 1);
    for i = 1 : L
        for j = 1 : L
            ind = (i-1) * L + j;
            ORG_STRUC.weight(ind) = (numIons(i) * numIons(j));
            S = S + (numIons(i) * numIons(j));
        end
    end
    
    ORG_STRUC.weight = ORG_STRUC.weight / S;
    
    % ionChange is just a handy little info thing
    for ind = 1:length(numIons)
        ORG_STRUC.ionChange(ind) = sum(numIons(1 : ind));
    end
    
    ORG_STRUC.ionCh = zeros(1, length(ORG_STRUC.ionChange) + 1);
    for i = 2 : length(ORG_STRUC.ionChange) + 1
        ORG_STRUC.ionCh(i) = ORG_STRUC.ionChange(i-1);
    end
end
%%%%%%%%%%%%%%%%%%%% END OF FIRST manipulations of the INPUT %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeORG()
%checkORG()
%[nothing,homePath] = unix('pwd');
%homePath(end) = [];
%checkPOTCARS(ORG_STRUC.atomType, ORG_STRUC.abinitioCode, ORG_STRUC.var_comp, homePath)
