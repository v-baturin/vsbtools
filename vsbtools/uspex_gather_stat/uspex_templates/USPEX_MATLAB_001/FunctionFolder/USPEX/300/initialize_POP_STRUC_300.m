function initialize_POP_STRUC_300()

% USPEX Version 9.4.0
% modulization

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, ...
'resFolder', {},'generation',{}, 'DoneOrder',{}, 'bodyCount', {}, 'ranking',{}, 'paretoFront',{},...
'bad_rank',{},  'CalcFold',{}, 'CalcFold_max',{}, 'fitness', {}, 'convex_hull',{}, 'finalOptimization', {});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'numIons', {}, 'INIT_numIons', {},...
    'struc_entr',{}, 'order',{},'dielectric_tensor',{}, 'gap',{}, 'hardness',{}, 'TE_property', {}, 'FINGERPRINT',{}, 'symg',{}, 'K_POINTS',{},...
    'mag_moment',{}, 'magmom_ions',{}, 'magmom_ini',{}, 'ldaU', {}, 'birefringence', {},...
    'S_order',{}, 'Step', {}, 'Enthalpies', {},'Fphon',{}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},'howCome',{},'JobID',{},'Folder',{}, 'Number',{});

POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;

%% create the composition file for 300
createCompostion_SingleBlock();

%create good initial population. Every individual fulfills hard constraints.
nSplit = length( ORG_STRUC.firstGeneSplit );
if nSplit == 0
    error('ERROR:  The minAt and maxAt is too small, please check your input file ...');
end

% The variable to store number of structures generated with wrong space group:
ORG_STRUC.wrong_spacegroups = 0;

% The variable defining number of structures in first generation to be generated with symmetry random
% the rest of structures wil be generated with topological random
howManyRandFirstGen = ORG_STRUC.fracRand/(ORG_STRUC.fracRand+ORG_STRUC.fracRandTop)*ORG_STRUC.initialPopSize;

for i = 1: ORG_STRUC.initialPopSize
    numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );
    numIons = ORG_STRUC.numIons*numBlocks;
    numMols = ORG_STRUC.numMols*numBlocks;
    if ~isempty(ORG_STRUC.numMols)  %% fragment
        %[candidate, lat] = Random_Init_310(i, ORG_STRUC.numMols);
        [candidate, lat] = Random_Init_310(i, numMols);
        POP_STRUC.POPULATION(i).howCome = '  Random  ';
    else
        %[candidate, lat] = Random_Init_300(i, ORG_STRUC.numIons);
        if i <= howManyRandFirstGen
            [candidate, lat] = Random_Init_300(i, numIons,0);%0 - symmetry based random
            POP_STRUC.POPULATION(i).howCome = '  Random  ';
        else
            [candidate, lat] = Random_Init_300(i, numIons,1);%1 - topological random
            POP_STRUC.POPULATION(i).howCome = '  RandTop  ';
        end
    end
    POP_STRUC.POPULATION(i).LATTICE = lat;
    POP_STRUC.POPULATION(i).COORDINATES = candidate;
    POP_STRUC.POPULATION(i).numIons = numIons;
end

if ORG_STRUC.wrong_spacegroups > 0
    disp([' ']);
    disp(['ATTENTION! In ' num2str(ORG_STRUC.wrong_spacegroups) ' / ' ...
        num2str(ORG_STRUC.initialPopSize) ' cases actually generated symmetry was different.']);
    disp([' ']);
end

if ORG_STRUC.spin == 1
    for i = 1: ORG_STRUC.initialPopSize
        POP_STRUC.POPULATION(i) = individual_Spin_Init(  POP_STRUC.POPULATION(i) );
    end
end

%%%%%%%%%%%% END defining the first generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% SEEDING %%%%%%%%%%%%%%%%%
pickUpSeeds();
if ORG_STRUC.doFing
    pickAntiSeeds();
end

Start_POP();


%% This inline function is created to

function createCompostion_SingleBlock()

global ORG_STRUC


if isempty(ORG_STRUC.minAt) || isempty(ORG_STRUC.maxAt)
    ORG_STRUC.minAt = sum(ORG_STRUC.numIons);
    ORG_STRUC.maxAt = sum(ORG_STRUC.numIons);
end

N_T = size(ORG_STRUC.numIons,1);
splitting = zeros(1,N_T);
findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
IPS = ORG_STRUC.initialPopSize;

fp = fopen('Seeds/compositions', 'w');
if exist('Seeds/Anti-compositions')
    movefile('Seeds/Anti-compositions', 'Seeds/Anti-compositions-back');
end
for i=1:size(ORG_STRUC.firstGeneSplit,1)
    for j=1:size(ORG_STRUC.firstGeneSplit,2)
        fprintf(fp, '%4d', ORG_STRUC.firstGeneSplit(i,j));
    end
    fprintf(fp, '\n');
end
fclose(fp);
