function MINHOP_Initialize()

%added USPEX_STRUC, lastly updated by Qiang Zhu (2014/02/18)
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

MINHOP_Initialize_POP();

POP_STRUC.resFolder = ORG_STRUC.resFolder;
%NEW We add USPEX_STRUC for analysis
USPEX_STRUC = struct('POPULATION',{}, 'GENERATION',{}, 'atomType', {}, 'Fp_weight', {});
USPEX_STRUC(1).GENERATION = struct('Best_enth', {}, 'Best_enth_relaxed',{}, 'BestID',{});
USPEX_STRUC.GENERATION(1) = QuickStart(USPEX_STRUC.GENERATION);
USPEX_STRUC(1).POPULATION = struct('COORDINATES', {}, 'LATTICE', {}, 'numIons',{}, ...
'symg',{}, 'howCome',{},'order',{}, 'FINGERPRINT', {}, 'K_POINTS', {}, ...
'ToCount',{}, 'gen',{}, 'Enthalpies', {}, 'Parents',{},'Vol',{},...
'PressureTensor',{}, 'coords0', {}, 'lat0',{}, 'PressureTensor0',{}, 'last_mut',{});

USPEX_STRUC.POPULATION(1) = QuickStart(USPEX_STRUC.POPULATION);
USPEX_STRUC.atomType = ORG_STRUC.atomType;
USPEX_STRUC.Fp_weight = ORG_STRUC.weight;

safesave('Current_POP.mat', POP_STRUC);
safesave([POP_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
POP_STRUC.resFolder = ORG_STRUC.resFolder;
if ORG_STRUC.varcomp == 1
   unixCmd(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals']);
   unixCmd(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals']);
   if ORG_STRUC.FullRelax > 0
   unixCmd(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals_relaxed']);
   unixCmd(['echo "Gen   ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals_relaxed']);
   end
else
   unixCmd(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals']);
   unixCmd(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals']);
   if ORG_STRUC.FullRelax > 0
   unixCmd(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/Individuals_relaxed']);
   unixCmd(['echo "Gen   ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY" >>' ORG_STRUC.resFolder '/BESTIndividuals_relaxed']);
   end
end

WriteGenerationStart();
