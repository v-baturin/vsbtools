function POP_STRUC  = NEB_creatPOPStruct()


global  ORG_STRUC

%---------------------------------------------------------------------------------------
memStruct_POP={
    'POPULATION'; 'resFolder'; 'step'; 'DoneOrder'; 'bodyCount'; 'genDone'; 'generation';
    };
memStruct_POPULATION={
    'Done'; 'ToDo'; 'JobID'; 'Folder'; 'Step';
    'LATTICE'; 'EPSILON';
    'CARTECOORDS'; 'COORDINATES'; 'caCoordVector'; 'fcCoordVector';
    'numIons';                              %USPEX compatiable
    'Energy'; 'Enthalpy'; 'VOLUME'; 'Press'; 'Stress'; 'Bulk'; 'preEnthalpy'; 'prePress'; 'preVOLUME'; 'prePathLength'
    'atomForcesMatrix'; 'cellStressMatirx'; %--- Read from the FILES
    'aCaForceMatrix'; 'aFrForceMatrix'; 'cCaStressMatrix'; 'cFrStressMatrix';
    'X_move'; 'u_move'; 'LATTICE_move'; 'EPSILON_move'; 'COORDINATES_move'; 'CARTECOORDS_move';
    'F_neb'; 'F_pro'; 'F_ela';
    'erraFelaRms'; 'erraFnebRms'; 'erraFproRms'; 'errcFnebRms'; 'errcFproRms'; 'errcFelaRms';
    'erraFelaMax'; 'erraFnebMax'; 'erraFproMax'; 'errcFnebMax'; 'errcFproMax'; 'errcFelaMax';
    'totFrForceVector'; 'totCaForceVector';
    'springK'; 'tuo'; 'pathLength'; 'freezing';
    'spacegroupNumber'; 'spacegroupSymbol';
    'Fire';
    'KPOINTS';
    'Error';
    'FINGERPRINT';'order';'struc_entr';'S_order';
    'symg';'Number';  %USPEX compatiable
    };

POP_STRUC = [];
for  i = 1:length(memStruct_POP)
    POP_STRUC = setfield(POP_STRUC, memStruct_POP{i}, {} );
end

POPULATION = [];
for  i = 1:length(memStruct_POPULATION)
    POPULATION = setfield(POPULATION, memStruct_POPULATION{i}, {} );
end
POP_STRUC(1).POPULATION = POPULATION;


POP_STRUC(1).POPULATION(1).Fire = struct ('dt', {}, 'a', {}, 'P', {}, 'v', {}, 'Ncount', {});

%-----------------------------------------------------------------

numImages = ORG_STRUC.numImages;
numIons = ORG_STRUC.numIons;
dimension = ORG_STRUC.dimension;
numDimension = 3*( sum(numIons) + dimension );

%------------
Fire = struct ('dt', 0, 'a', 0, 'P', 0, 'v', zeros(numDimension,1), 'Ncount', 0);

%   'POPULATION'; 'resFolder'; 'step'; 'DoneOrder'; 'bodyCount'; 'genDone'; 'Done'; 'ToDo'; 'JobID'; 'Folder';
POP_STRUC.resFolder = [];
POP_STRUC.genDone   = 0;
POP_STRUC.step      = 0;
POP_STRUC.generation= POP_STRUC.step;
POP_STRUC.bodyCount = 0;
POP_STRUC.DoneOrder = zeros(numImages,1);

for i = 1:numImages
    POP_STRUC.POPULATION(i).Folder = 0;
    POP_STRUC.POPULATION(i).Done = 0;
    POP_STRUC.POPULATION(i).ToDo = 1;
    %   .LATTICE'; .EPSILON';
    POP_STRUC.POPULATION(i).LATTICE = zeros(3,3);
    POP_STRUC.POPULATION(i).EPSILON = zeros(3,3);
    %  'numIons', 'Energy'; 'Enthalpy'; 'VOLUME'; 'Press'; 'Stress'; 'Bulk'; 'prePress'; 'preVOLUME';
    POP_STRUC.POPULATION(i).numIons= ORG_STRUC.numIons;
    POP_STRUC.POPULATION(i).Energy = 0;
    POP_STRUC.POPULATION(i).Enthalpy = 0;
    POP_STRUC.POPULATION(i).VOLUME   = 0;
    POP_STRUC.POPULATION(i).Press    = 0;
    POP_STRUC.POPULATION(i).Stress   = zeros(3,3);
    POP_STRUC.POPULATION(i).Bulk     = 0;
    POP_STRUC.POPULATION(i).prePress = 0;
    POP_STRUC.POPULATION(i).preVOLUME = 0;
    %  'atomForcesMatrix'; 'cellStressMatirx';
    POP_STRUC.POPULATION(i).atomForcesMatrix = zeros(sum(numIons),3);
    POP_STRUC.POPULATION(i).cellStressMatirx = zeros(3,3);
    % 'CARTECOORDS'; 'COORDINATES'; 'caCoordVector'; 'fcCoordVector';
    POP_STRUC.POPULATION(i).CARTECOORDS= zeros(sum(numIons),3);
    POP_STRUC.POPULATION(i).COORDINATES= zeros(sum(numIons),3);
    POP_STRUC.POPULATION(i).caCoordVector = zeros( numDimension, 1 );
    POP_STRUC.POPULATION(i).frCoordVector = zeros( numDimension, 1 );
    
    %  'totFrForceVector'; 'totCaForceVector';
    POP_STRUC.POPULATION(i).totCaForceVector = zeros( numDimension, 1 );
    POP_STRUC.POPULATION(i).totFrForceVector = zeros( numDimension, 1 );
    
    % 'springK'; 'tuo'; 'pathLength'; 'freezing';
    POP_STRUC.POPULATION(i).tuo       = zeros(numDimension,1);
    POP_STRUC.POPULATION(i).springK   = 0;
    POP_STRUC.POPULATION(i).Bulk      = 0;
    POP_STRUC.POPULATION(i).freezing  = 0;
    POP_STRUC.POPULATION(i).pathLength= 0;
    % 'F_neb'; 'F_pro'; 'F_ela';
    POP_STRUC.POPULATION(i).F_neb = zeros( numDimension, 1 );
    POP_STRUC.POPULATION(i).F_pro = zeros( numDimension, 1 );
    POP_STRUC.POPULATION(i).F_ela = zeros( numDimension, 1 );
    %  'erraFelaRms';  'erraFnebRms';  'erraFproRms';  'errcFnebRms';  'errcFproRms';  'errcFelaRms';
    %  'erraFelaMax'; 'erraFnebMax'; 'erraFproMax'; 'errcFnebMax'; 'errcFproMax'; 'errcFelaMax';
    POP_STRUC.POPULATION(i).errcFnebMax = 0;
    POP_STRUC.POPULATION(i).errcFproMax = 0;
    POP_STRUC.POPULATION(i).errcFelaMax = 0;
    POP_STRUC.POPULATION(i).erraFnebMax = 0;
    POP_STRUC.POPULATION(i).erraFproMax = 0;
    POP_STRUC.POPULATION(i).erraFelaMax = 0;
    POP_STRUC.POPULATION(i).errcFnebRms = 0;
    POP_STRUC.POPULATION(i).errcFproRms = 0;
    POP_STRUC.POPULATION(i).errcFelaRms = 0;
    POP_STRUC.POPULATION(i).erraFnebRms = 0;
    POP_STRUC.POPULATION(i).erraFproRms = 0;
    POP_STRUC.POPULATION(i).erraFelaRms = 0;
    
    %    POP_STRUC.POPULATION(i).Hessian = zeros(numDimension, numDimension);
    %   'X_move'; 'u_move'; .LATTICE_move'; .EPSILON_move'; 'COORDINATES_move'; 'CARTECOORDS_move';
    POP_STRUC.POPULATION(i).X_move = zeros(numDimension,1);
    POP_STRUC.POPULATION(i).u_move = zeros(numDimension,1);
    POP_STRUC.POPULATION(i).LATTICE_move = zeros(3,3);
    POP_STRUC.POPULATION(i).EPSILON_move = zeros(3, 3);
    POP_STRUC.POPULATION(i).CARTECOORDS_move = zeros(sum(numIons),3);
    POP_STRUC.POPULATION(i).COORDINATES_move = zeros(sum(numIons),3);
    %   'spacegroupNumber'; 'spacegroupSymbol';
    POP_STRUC.POPULATION(i).spacegroupNumber  = 0;
    POP_STRUC.POPULATION(i).spacegroupSymbol = [];
    POP_STRUC.POPULATION(i).symg             = 0;
    %   'KPOINTS';
    POP_STRUC.POPULATION(i).KPOINTS = zeros(1,3);
    %  'Fire'
    POP_STRUC.POPULATION(i).Fire = Fire;
    POP_STRUC.POPULATION(i).Step = 1;
    POP_STRUC.POPULATION(i).Error = 0;
    % USPEX compatiable
    POP_STRUC.POPULATION(i).Number = 0;
end
