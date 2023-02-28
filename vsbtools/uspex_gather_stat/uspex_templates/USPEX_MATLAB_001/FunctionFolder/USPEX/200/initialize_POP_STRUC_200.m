function initialize_POP_STRUC_200()
%This routine is to initiallize surface - 200
%Last updated by Qiang Zhu, 2013/10/09

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, ...
'resFolder', {},'generation',{}, 'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},...
'bad_rank',{},  'CalcFold',{}, 'CalcFold_max',{}, 'fitness', {}, 'convex_hull',{});


POP_STRUC(1).POPULATION = struct('COORDINATES',{}, 'LATTICE',{},'numIons',{},'chanAList',{},'cell',{},...
                         'INIT_COORD',{},'INIT_LAT',{},'INIT_numIons',{},'Enthalpies',{},'S_order',{},...
                            'struc_entr', {},'order',{}, 'FINGERPRINT', {}, 'K_POINTS', {},'Step', {},...
                     'Error',{},'Done',{}, 'ToDo',{},'Parents',{},'howCome',{},'JobID',{},'Folder',{},...
              'Surface_COORDINATES',{}, 'Surface_LATTICE',{}, 'Surface_order',{},'Surface_numIons',{},...
                    'Bulk_COORDINATES',{},'Bulk_LATTICE',{}, 'Bulk_numIons',{},'Number',{}, 'symg',{},...
                                       'mag_moment',{}, 'magmom_ions',{}, 'magmom_ini',{}, 'ldaU', {});

POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{});


POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;

Vacuum   = ORG_STRUC.vacuumSize(1);
CellList = findcell(ORG_STRUC.bulk_lat,ORG_STRUC.reconstruct);
fp = fopen('Seeds/cell', 'w');
for i=1:size(CellList,1)
    fprintf(fp, '%4d %4d %4d %4d\n', CellList(i,:));
end
fclose(fp);
   
for goodPop = 1: ORG_STRUC.initialPopSize
    ID = ceil(rand()*size(CellList,1));
    cell=CellList(ID,:);
    ANS  = Random_Init_200(goodPop, Vacuum, cell);
    POP_STRUC.POPULATION(goodPop).COORDINATES         = ANS(1).candidate;
    POP_STRUC.POPULATION(goodPop).numIons             = ANS(1).numIons;
    POP_STRUC.POPULATION(goodPop).LATTICE             = ANS(1).lat;
    POP_STRUC.POPULATION(goodPop).chanAList           = ANS(1).chanAList;
    POP_STRUC.POPULATION(goodPop).Surface_LATTICE     = ANS(1).sur_lat;
    POP_STRUC.POPULATION(goodPop).Surface_COORDINATES = ANS(1).sur_candidate;
    POP_STRUC.POPULATION(goodPop).Surface_numIons     = ANS(1).sur_numIons;
    POP_STRUC.POPULATION(goodPop).Bulk_LATTICE        = ANS(1).bulk_lat;
    POP_STRUC.POPULATION(goodPop).Bulk_COORDINATES    = ANS(1).bulk_pos;
    POP_STRUC.POPULATION(goodPop).Bulk_numIons        = ANS(1).bulk_numIons;
    POP_STRUC.POPULATION(goodPop).cell                = cell;
    POP_STRUC.POPULATION(goodPop).howCome             = '  Random  ';
end

if ORG_STRUC.spin == 1
    for i = 1: ORG_STRUC.initialPopSize
        POP_STRUC.POPULATION(i) = individual_Spin_Init(  POP_STRUC.POPULATION(i) );
    end
end

pick_Seeds_surface();
Start_POP();
