function META_Initialize_POP()

global ORG_STRUC
global POP_STRUC

POP_STRUC = struct('POPULATION',{}, 'resFolder', {},'generation',{}, 'DoneOrder',{}, 'basicStructureNumber',{},...
     'bodyCount', {}, 'finalOptimization', {}, 'ranking',{}, 'lat',{}, 'lat_6',{}, 'usedSoftModes', {});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'FINGERPRINT', {}, 'K_POINTS', {},'Step', {}, 'Enthalpies', {},...
    'coords0', {}, 'lat0', {}, 'PressureTensor0', {},...    % parameters before 'final' optimization
    'Done',{},'ToDo',{}, 'JobID',{},'Folder',{}, 'numIons', {}, 'eignFre',{}, 'eignVec',{}, 'eignSupercells',{}, 'Number',{},'Error',{},...
    'Softmode_Fre', {}, 'Softmode_num',{}, 'PressureTensor', {}, 'Parents', {}, 'order',{}, 'superCell',{}, 'symg',{});

POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC.generation = 0;
POP_STRUC.bodyCount = 0;
POP_STRUC.usedSoftModes = [];
POP_STRUC.lat   = ORG_STRUC.lattice;
POP_STRUC.lat_6 = MattoVec(ORG_STRUC.lattice);

POP_STRUC.POPULATION(1).superCell = [1 1 1]; % describes the supercell regarding the initial structure
POP_STRUC.POPULATION(1).COORDINATES = ORG_STRUC.coordinates;
POP_STRUC.POPULATION(1).numIons = ORG_STRUC.numIons;
disp('Generation 0 is initialized.');
disp(' ');
META_Start_POP();
