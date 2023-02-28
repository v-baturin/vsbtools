function initialize_POP_STRUC_000()

% USPEX Version 9.4.0
% modulization

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, ...
'resFolder', {},'generation',{}, 'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},...
'bad_rank',{},  'CalcFold',{}, 'CalcFold_max',{}, 'fitness', {}, 'convex_hull',{});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {}, 'numIons',{}, ...
'Vol',{},'struc_entr',{},'order',{},'FINGERPRINT',{},'K_POINTS',{},'Step',{}, 'Enthalpies',{}, 'Error',{},'Done',{},'ToDo',{},'Parents',{},...
'mag_moment',{}, 'magmom_ions',{}, 'magmom_ini',{}, 'ldaU', {}, ...
'S_order',{}, 'howCome',{},'JobID',{},'Folder',{}, 'Number',{}, 'symg', {});
POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;

numIons = ORG_STRUC.numIons;

for i = 1: ORG_STRUC.initialPopSize
    [candidate, lat] = Random_Init_000(i, numIons);
    POP_STRUC.POPULATION(i).LATTICE = lat;
    POP_STRUC.POPULATION(i).COORDINATES = candidate;
    POP_STRUC.POPULATION(i).numIons = numIons;
    POP_STRUC.POPULATION(i).howCome = '  Random  ';;
end


if ORG_STRUC.spin == 1
    for i = 1: ORG_STRUC.initialPopSize
        POP_STRUC.POPULATION(i) = individual_Spin_Init(  POP_STRUC.POPULATION(i) );
    end
end

%%%%%%%%%%%%%%%%%%%% SEEDING %%%%%%%%%%%%%%%%%
pickUpSeeds();
if ORG_STRUC.doFing
 pickAntiSeeds();
end
Start_POP();
