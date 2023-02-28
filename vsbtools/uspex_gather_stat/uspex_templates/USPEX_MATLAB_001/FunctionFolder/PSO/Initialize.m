function Initialize()

%added USPEX_STRUC, lastly updated by Qiang Zhu (2014/02/18)
global ORG_STRUC
global POP_STRUC
global PSO_STRUC


if ~exist('Seeds','dir')
    unixCmd('mkdir Seeds');
end
if ~exist('AntiSeeds','dir')
    unixCmd('mkdir AntiSeeds');
end

Initialize_POP_STRUC_PSO();

POP_STRUC.resFolder = ORG_STRUC.resFolder;
%NEW We add USPEX_STRUC for analysis
PSO_STRUC = struct('POPULATION',{}, 'GENERATION',{});
PSO_STRUC(1).GENERATION = struct('quasiEntropy', {}, 'convex_hull', {}, 'composEntropy',{}, 'BestID',{});
PSO_STRUC.GENERATION(1) = QuickStart(PSO_STRUC.GENERATION);
PSO_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons',{}, ...
    'symg',{}, 'howCome',{},'order',{}, 'Fitness', {}, 'FINGERPRINT', {}, 'K_POINTS', {}, ...
    'ToCount',{},'S_order',{}, 'gen',{}, 'struc_entr',{}, 'Enthalpies', {}, 'Parents',{},'Vol',{},...
    'cell',{},'Surface_numIons',{},'numBlocks',{},'gap','hardness','mag_moment','dielectric_tensor');
PSO_STRUC.POPULATION(1) = QuickStart(PSO_STRUC.POPULATION);

safesave('Current_POP.mat', POP_STRUC);
safesave([POP_STRUC.resFolder '/PSO.mat'], PSO_STRUC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unixCmd(['echo "Gen   ID    Origin   Composition    Enthalpy   Volume  Density   Fitness   KPOINTS  SYMM  Q_entr A_order S_order" >>' ORG_STRUC.resFolder '/Individuals']);
unixCmd(['echo "                                      (eV)     (A^3)   (g/cm^3)" >>' ORG_STRUC.resFolder '/Individuals']);

unixCmd(['echo "  ID   Origin     Composition  Enthalpy(eV)  Volume(A^3)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/OUTPUT.txt']);
unixCmd(['echo " ID    Origin    Enthalpy   Parent-E   Parent-ID" >>' ORG_STRUC.resFolder '/origin']);
WriteGenerationStart();
