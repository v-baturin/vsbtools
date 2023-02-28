function Initialize()

%added USPEX_STRUC, lastly updated by Qiang Zhu (2014/02/18)
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
global POOL_STRUC


[homePath, USPEXPath] = workingPath();

if ~exist('Seeds','dir')
    unixCmd('mkdir Seeds');
end
if ~exist('AntiSeeds','dir')
    unixCmd('mkdir AntiSeeds');
end


if ORG_STRUC.fixRndSeed > 0
	rng( ORG_STRUC.fixRndSeed+1, 'twister' );
end

% Preprocess data before initialization of a particular calculation type:
if (ORG_STRUC.dimension==3)  &&  (ORG_STRUC.varcomp==1)
    N_T = size(ORG_STRUC.numIons,1);
    splitting = zeros(1,N_T);
    findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
end

% Initialization part:
calcType_str = [num2str(ORG_STRUC.dimension) '' num2str(ORG_STRUC.molecule) '' num2str(ORG_STRUC.varcomp)];
if str2num(calcType_str) < 0
    calcType = -1*str2num(calcType_str);
    calcType_str = ['M' num2str(calcType)];
end

%NEW We add USPEX_STRUC for analysis
USPEX_STRUC = struct('POPULATION',{}, 'GENERATION',{}, 'SYSTEM', {});

USPEX_STRUC(1).GENERATION = struct('quasiEntropy', {}, 'convex_hull', {}, 'composEntropy',{}, 'BestID',{}, 'ID',{});
USPEX_STRUC.GENERATION(1) = QuickStart(USPEX_STRUC.GENERATION);
USPEX_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons',{}, 'MOLECULES',{},...
'symg',{}, 'howCome',{},'order',{}, 'Fitness', {}, 'FINGERPRINT', {}, 'K_POINTS', {}, ...
'ToCount',{},'S_order',{}, 'gen',{}, 'struc_entr',{}, 'Enthalpies', {}, 'Parents',{},'Vol',{},...
'cell',{},'Surface_numIons',{},'numBlocks',{},'Dis_ConvexHull',{},'gap',{},'Fphon',{},'hardness',{},...
'mag_moment',{},'magmom_ions',{}, 'magmom_ini',{}, 'ldaU',{}, 'dielectric_tensor',{}, 'birefringence',{});
USPEX_STRUC(1).SYSTEM    = struct( 'atomType', {}, 'atomType_symbol',{}, 'Fp_weight', {}, 'STDMOL',{});

USPEX_STRUC.POPULATION(1) = QuickStart(USPEX_STRUC.POPULATION);
%%%%%%%%%%%%%%%%% Initialize the global parameter for data analysis
atomType = ORG_STRUC.atomType;
USPEX_STRUC.SYSTEM(1).atomType = atomType;
for i = 1:length(atomType)
   USPEX_STRUC.SYSTEM(1).atomType_symbol{i} = megaDoof(atomType(i));
end
USPEX_STRUC.SYSTEM(1).Fp_weight = ORG_STRUC.weight;
USPEX_STRUC.SYSTEM(1).STDMOL    = ORG_STRUC.STDMOL;


%Extention for topology based random structure generation
USPEX_STRUC(1).EXT_RandTop = struct('STRUCTURES',{},'COORDINATING_NUMBERS',{},'CN_grid',{});
USPEX_STRUC.EXT_RandTop(1).STRUCTURES = struct('LATTICE',{},'COORDINATES',{});
USPEX_STRUC.EXT_RandTop(1).COORDINATING_NUMBERS = struct('COMPOSITION',{},'CN',{});
USPEX_STRUC.EXT_RandTop(1).CN_grid = cell(0);

if ORG_STRUC.dimension==3
    if isempty(ORG_STRUC.maxAt)
        grid_sz = sum(ORG_STRUC.numIons)+1;
    else
        grid_sz = ORG_STRUC.maxAt+1;
    end

    for i = 1:size(atomType,2)
        a_size = repmat([grid_sz],1,size(atomType,2));
        a_grid = ceil(repmat([3],a_size));
        USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = a_grid;
        a_grid = ceil(repmat([12],a_size));
        USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = a_grid;
    end
end


path(path,[USPEXPath '/FunctionFolder/USPEX/' calcType_str]);
eval(['initialize_POP_STRUC_' calcType_str '()']);
POP_STRUC.resFolder = ORG_STRUC.resFolder;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

safesave('Current_POP.mat', POP_STRUC);
safesave([POP_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
if exist('POOL_STRUC','var')
    safesave([POP_STRUC.resFolder '/POOL.mat'], POOL_STRUC);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unixCmd(['echo "Gen   ID    Origin   Composition    Enthalpy   Volume  Density   Fitness   KPOINTS  SYMM  Q_entr A_order S_order Magmoment-Type" >>' ORG_STRUC.resFolder '/Individuals']);
if (ORG_STRUC.dimension==-4)  &&  (ORG_STRUC.molecule==0)  &&  (ORG_STRUC.varcomp==0)
   [a,b]= unix(['echo "                                   (kcal/mol)  (A^3)   (g/cm^3)" >>' ORG_STRUC.resFolder '/Individuals']);
   [a,b]= unix(['echo "  ID   Origin     Composition  Enthalpy      Volume(A^3)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
   [a,b]= unix(['echo "                              (kcal/mol)" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
else
   [a,b]= unix(['echo "                                      (eV)     (A^3)   (g/cm^3)" >>' ORG_STRUC.resFolder '/Individuals']);
   [a,b]= unix(['echo "  ID   Origin     Composition  Enthalpy(eV)  Volume(A^3)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
end

[a,b]=unix(['echo " ID    Origin    Enthalpy   Parent-E   Parent-ID" >>' ORG_STRUC.resFolder '/origin']);
WriteGenerationStart();
