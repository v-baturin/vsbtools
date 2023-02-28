function initialize_POP_STRUC_M200()
% 2D

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, ...
'resFolder', {},'generation',{}, 'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},...
'bad_rank',{},  'CalcFold',{}, 'CalcFold_max',{}, 'fitness', {}, 'convex_hull',{});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {}, 'numIons',{},...
'struc_entr',{}, 'order',{},'FINGERPRINT',{},'K_POINTS',{},'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},...
'S_order',{},'howCome',{},'JobID',{},'Folder',{},'COORDINATES_2D',{},'LATTICE_2D',{},'Number',{}, 'symg',{},...
'mag_moment',{}, 'magmom_ions',{}, 'magmom_ini',{}, 'ldaU', {});


POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},...
                         'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;

% create the composition file for 300
nSplit = createCompostion_SingleBlock();
if nSplit == 0
   disp('ERROR:  The minAt and maxAt is too small, please check your input file ...');
end

%create good initial population. Every individual fulfills hard constraints.
for i = 1: ORG_STRUC.initialPopSize
   numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand*nSplit) );
    if ~isempty(ORG_STRUC.numMols)
      numMols = ORG_STRUC.numMols*numBlocks;
      [candidate_2D, lat_2D, candidate, lat] = Random_Init_M210(i, numMols);
      numIons = ORG_STRUC.numIons;
    else
      numIons = ORG_STRUC.numIons*numBlocks;
      [candidate_2D, lat_2D, candidate, lat] = Random_Init_M200(i, numIons);
    end
      POP_STRUC.POPULATION(i).LATTICE = lat;
      POP_STRUC.POPULATION(i).LATTICE_2D = lat_2D;
      POP_STRUC.POPULATION(i).COORDINATES = candidate;
      POP_STRUC.POPULATION(i).COORDINATES_2D = candidate_2D;
      POP_STRUC.POPULATION(i).howCome = '  Random  ';
      POP_STRUC.POPULATION(i).numIons = numIons;
end

if ORG_STRUC.spin == 1
    for i = 1: ORG_STRUC.initialPopSize
        POP_STRUC.POPULATION(i) = individual_Spin_Init(  POP_STRUC.POPULATION(i) );
    end
end

    pick_Seeds();
Start_POP();

%% This inline function is created to

function nSplit = createCompostion_SingleBlock()

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
nSplit = size(ORG_STRUC.firstGeneSplit,1);
