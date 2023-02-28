function initialize_POP_STRUC_M300()

global ORG_STRUC
global POP_STRUC
POP_STRUC= struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{},...
        'resFolder', {},'generation',{}, 'DoneOrder',{}, 'bodyCount', {}, ...
          'ranking',{}, 'bad_rank',{},  'CalcFold',{}, 'CalcFold_max',{},...
                                         'fitness', {}, 'convex_hull',{});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, ...
    'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {}, 'numIons', {}, ...
      'struc_entr',{}, 'order',{}, 'FINGERPRINT', {}, 'K_POINTS', {},...
        'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},...
     'Parents',{}, 'S_order',{}, 'howCome',{},'JobID',{},'Folder',{},...
                'typesAList',{}, 'chanAList',{}, 'GB_COORDINATES',{}, ...
   'GB_LATTICE',{}, 'GB_order',{},'GB_numIons',{},'GB_typesAList',{},...
          'Bulk_COORDINATES',{},'Bulk_LATTICE',{}, 'Bulk_numIons',{},...
                        'Bulk_typesAList',{}, 'Number',{}, 'symg',{});

POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{}, ...
                       'fingerprint',{},'eignFre',{},'eignVec',{}, ...
                              'Softmode_Fre',{},'Softmode_num',{});


POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;
%% create the composition file for 300
createCompostion_SingleBlock();
nSplit = length( ORG_STRUC.firstGeneSplit );
if nSplit == 0
    error('ERROR: minAt and maxAt is too small, check your input file ...');
end

for i = 1:ORG_STRUC.initialPopSize
    numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );
    numIons = ORG_STRUC.numIons*numBlocks;

    POP = Random_Init_M300(i, numIons);
    POP_STRUC.POPULATION(i).Bulk_LATTICE     = POP.Bulk_LATTICE;
    POP_STRUC.POPULATION(i).Bulk_COORDINATES = POP.Bulk_COORDINATES;
    POP_STRUC.POPULATION(i).Bulk_typesAList  = POP.Bulk_typesAList;
    POP_STRUC.POPULATION(i).Bulk_numIons     = POP.Bulk_numIons;
    POP_STRUC.POPULATION(i).COORDINATES      = POP.COORDINATES;
    POP_STRUC.POPULATION(i).numIons          = POP.numIons;         
    POP_STRUC.POPULATION(i).LATTICE          = POP.LATTICE; 
    POP_STRUC.POPULATION(i).typesAList       = POP.typesAList;
    POP_STRUC.POPULATION(i).GB_LATTICE       = POP.GB_LATTICE;      
    POP_STRUC.POPULATION(i).GB_COORDINATES   = POP.GB_COORDINATES;
    POP_STRUC.POPULATION(i).GB_numIons       = POP.GB_numIons;
    POP_STRUC.POPULATION(i).GB_typesAList    = POP.GB_typesAList;  
    POP_STRUC.POPULATION(i).chanAList        = POP.chanAList;
    POP_STRUC.POPULATION(i).howCome          = ' Random ';

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
