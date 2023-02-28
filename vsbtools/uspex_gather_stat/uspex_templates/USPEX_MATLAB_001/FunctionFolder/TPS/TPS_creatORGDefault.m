function memStruct = TPS_creatORGDefault() 


global ORG_STRUC



%---------------------------------------------
memStruct.numeral = {
    'orderParaType';
    'shiftRatio';
    'fixRandSeeds';
    'maxMDLoop';
    'numIterations';
    'SGtolerance';
    'continuedNow';
    'pickUpYN';     'pickUpGen'; 'pickUpFolder'; 
    'submitCount';  
    'remoteRegime';              'numProcessors';     'numParallelCalcs';  
    'UnitType';  'Ry2eV';        'KBar2GPa';  'fs';            'bohr2A';        'amu2eVs2_A2';   'THz2cm_1'; 
    'read_Handles';  'write_Handles'; 
    'platform';
    'repeatForStatistics';
    'maxErrors';
    };

memStruct.array = {
    'numSpecies';  'numIons';  'mass';  'opCriteria';
    'atomType';     'atomSymbol'; 'speciesSymbol';
    'atomicMasses'
    'amplitudeShoot';     'magnitudeShoot'; 
    'SuperCell';
    'Kresol'; 
    'abinitioCode';
    'ExternalPressure';
    };

memStruct.string = {
    'username';  
    'whichCluster';   'wallTime';
    'remoteFolder';   'specificFolder'; 
    'inputFile';
    'commandExecutable';
    'cmdCalcOp';
    'cmdCalcHT';
    'orderParameterFile'; 'enthalpyTemptureFile';
    'trajectoryFile';  'MDrestartFile';
    'resFolder';     'homePath'; 'USPEXPath';
    'log_file';
   };


ORG_STRUC = struct( 'dimension', 3 );


for  i = 1:length(memStruct.numeral)
    ORG_STRUC = setfield(ORG_STRUC, memStruct.numeral{i}, 0 );
end
for  i = 1:length(memStruct.array)
    ORG_STRUC = setfield(ORG_STRUC, memStruct.array{i}, '' );
end
for  i = 1:length(memStruct.string)
    ORG_STRUC = setfield(ORG_STRUC, memStruct.string{i}, '' );
end

%---------------------------------------------
%%%%%%%%% Create folders and stuff %%%%%%%%%%%%%%%%

[ORG_STRUC(1).homePath, ORG_STRUC(1).USPEXPath] = workingPath(); 

% ****************************** %
% ****************************** %
% *         Platform           * %
% *      0: Nonparallel        * %        
% *      1: local              * %    
% *      2: remote             * %
% *      3: CFN                * %
% *      4: QSH                * %
% *      5: QSH2               * % 
% *      6: xservDE            * %
% *      7: MIPT               * %    
ORG_STRUC.platform = 0;    


% *****************************************
% *****************************************
% *                System                 *
% *****************************************
% *****************************************


%-- ExternalPressure (GPa)   --------------
%   ExternalPressure = 0 
%------------------------------------------ 
ORG_STRUC.ExternalPressure = 0;

%-- Temperature (K)  ----------------------
%   Temperature = 0
%------------------------------------------
ORG_STRUC.Temperature = 0;

%-- SuperCell -----------------------------
%  SuperCell = [ 1 1 1 ]
%------------------------------------------
ORG_STRUC.SuperCell = [1 1 1];

ORG_STRUC.numIterations = 1000;
ORG_STRUC.maxMDLoop = 10;
ORG_STRUC.maxErrors = 5;

% *****************************************
% *****************************************
% *               Shooting                *
% *****************************************
% *****************************************

ORG_STRUC.amplitudeShoot = [  0.1, 0.1  ];
ORG_STRUC.magnitudeShoot = [ 1.05, 1.05 ];
ORG_STRUC.shiftRatio     = 0.1;


ORG_STRUC.orderParameterFile   = 'fp.dat';
ORG_STRUC.enthalpyTemptureFile = 'HT.dat';
ORG_STRUC.trajectoryFile       = 'traj.dat';
%ORG_STRUC.MDrestartFile        ='';


% *****************************************
% *****************************************
% *               OUTPUT                  *
% *****************************************
% *****************************************

%-- Structral OUTPUT Format ----------
%   1. xsf  xcrystal         ( default )
%   2. VASP POSCAR           
%   3. XYZ  fromat                   
%------------------------------------------

%-- Unit Type -----------------------------
%   1. A/eV/GPa/ps      for VASP  (default)
%   2. bohr/eV/kBar/ps  for QE
%------------------------------------------
ORG_STRUC.UnitType = 1;

%-- Pathway Print File Steps --------------
%   PrintStep = 1
%------------------------------------------
ORG_STRUC.PrintStep = 1; 

% *****************************************
% *****************************************
% *   DETAILS OF AB INITIO CALCULATIONS   * 
% *****************************************
% *****************************************

% ****************************** %
% ****************************** %
% *        Abinit code         * %
% ****************************** %
% 1:vasp 2:siesta 3:gulp 4:lammps 5:NN 6:dmacrys 7:CP2K, 8:QE 9:FHIaims: 10:ATK  11:CASTEP 
%
%-Current Supported Code in VCNEB
%    1: vasp, 3:gulp, 8:QE 
%-Further Plan:
%    7: CP2K
ORG_STRUC.abinitioCode = 1; 


%-- K Points Resolution ------------------
%   kresol = 0.075           (default)
%-----------------------------------------
ORG_STRUC.Kresol = 0.075;


% ****************************** % 
% ****************************** %
% *      Space Group           * %  
% ****************************** %
% ****************************** %
% Output symmetrized cif file
%ORG_STRUC.symmetrize = 0;
% Space group determination tolerance
ORG_STRUC.SGtolerance = 0.10;

ORG_STRUC.pickUpYN = 0;


%-- Some Constants:
ORG_STRUC.dimension  = 3;
ORG_STRUC.bohr2A     = 0.529177249;
ORG_STRUC.Ry2eV      = 13.6056981;
ORG_STRUC.KBar2GPa   = 10;
ORG_STRUC.fs         = 10^(-15);
ORG_STRUC.amu2eVs2_A2= 1.0364269*10^(-28);
ORG_STRUC.THz2cm_1   = (10^10)/29979245800 ;

%-----------------------------------------------

ORG_STRUC.log_file = 'OUTPUT.txt';
ORG_STRUC.specificFolder = 'Specific';

%% for compatibility 
ORG_STRUC.varcomp = 0;  
ORG_STRUC.molecule= 0; 
