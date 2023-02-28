function memStruct = NEB_creatORGDefault()


global ORG_STRUC



%---------------------------------------------
memStruct.numeral = {
    'numImages';
    'CalcType';      'optRelaxType';   'optVarImage';      'optVarK';
    'optimizerType';  'optReadImages'; 'optMethodCIDI';    'optFreezing';
    'ionCh';         'whetherConstraint';
    'dimension';
    'Temperature';   'SuperCell'; 'ExternalPressure';
    'dt';           'K_min';         'K_max';   'Kconstant';
    'VarPathLength';
    'startCIDIStep';
    'ConvThreshold';
    'numSteps';   'FormatType';      'PrintStep';
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
    'RmaxFing'; 'deltaFing'; 'sigmaFing'; 'toleranceFing';
    };

memStruct.array = {
    'atomType';     'atomSymbol';
    'numIons';      'numSpecies';
    'SuperCell';
    'whichCI';      'whichDI';     'pickupImages';
    'Kresol';
    'abinitioCode';
    };

memStruct.string = {
    'username';
    'whichCluster';   'wallTime';
    'remoteFolder';   'specificFolder';
    'inputFile';
    'commandExecutable';
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

% ****************************** %
% ****************************** %
% *    Calculation Option      * %
% ****************************** %
% ****************************** %

%--  Calculation Type -------------------
%    1: VC-NEB               (default)
%    2: Structrual Relax
%    3: Structural Relax + T
%----------------------------------------
ORG_STRUC.CalcType = 1;

%-- Relaxation Type ---------------------
%   1: Only Relax Atom
%   2: Only Relax Cell
%   3: Relax      Cell+Atom  (default)
%---------------------------------------
ORG_STRUC.optRelaxType = 3;

%-- Strucutre Optimization Algorithm  ---
%    1: Quick-MIN / QM
%    2: FIRE                      (default)
%    3: BFGS                      (not yet)
%    4: Fast Method for Cell Opt  (not yet)
%-----------------------------------------
ORG_STRUC.optimizerType = 1;

%-- Method for Reading the Images  -------
%   0: all Images Read
%   1: only the First&Last Images (default)
%   2: auotmatic interpolation
%------------------------------------------
ORG_STRUC.optReadImages = 1;

%-- Number of the Images ------------------
%   numImages  =  1
%------------------------------------------
ORG_STRUC.numImages = 1;

%-- Number of the Calculation Steps -------
%   numSteps   = 10
%------------------------------------------
ORG_STRUC.numSteps = 10;

%-- Convergence Tolerance -----------------
%  ConvThreshold = 1e-3;
%------------------------------------------
ORG_STRUC.ConvThreshold = 3e-3;

%-- Relaxation Time Step  -----------------
%  dt = 0.05;
%------------------------------------------
ORG_STRUC.dt = 0.05;

% ****************************** %
% ****************************** %
% *       VCNEB Option         * %
% ****************************** %
% ****************************** %

%-- Spring Constant Method  ---------------
%    1: Constant-K              (default)
%    2: Variable-K
%------------------------------------------
ORG_STRUC.optVarK = 1;


%-- Variable Image Number Method  ---------------
%    1: Constant-Image number
%    1: Variable-Image number   (default)
%------------------------------------------
ORG_STRUC.optVarImage = 1;
ORG_STRUC.VarPathLength = 0;

%-- Spring Constant Value   ---------------
%   K_min = 5;
%   K_max = 5;
%------------------------------------------
ORG_STRUC.K_min = 5;
ORG_STRUC.K_max = 5;

ORG_STRUC.Kconstant=5;

%-- Image Freezing Option  ----------------
%   0: NOT Freezing               (default)
%   1: Freezing after FreezingSteps
%------------------------------------------
ORG_STRUC.optFreezing = 0;

%-- CI,DI Method --------------------------
%    0: no CI,DI Method           (default)
%    1: singleCI Method (TS)
%   -1: singleDI Method (minimal search)
%    2: multiCI-DI Method         (not yet)
%-----------------------------------------
ORG_STRUC.optMethodCIDI = 0;

ORG_STRUC.whichCI = 0;
ORG_STRUC.whichDI = 0;

%ORG_STRUC.pickupImages = [];

%-- Climbing Steps for VCNEB Calculatio --
%   startCIDIStep = 100000
%------------------------------------------
ORG_STRUC.startCIDIStep = 100000;

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

%%--- Must Read Parameters !-----------------
% 1) numIons
% 2) atomType
%--------------------------------------------

% ORG_STRUC.startCIDIStep = str2num(startCIDIStep);
% ORG_STRUC.pickupImages = str2num(pickupImages);


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


%--- Energy tag ---%
ORG_STRUC.E_IniImage = 10E+9;
ORG_STRUC.fp_IniImage=[];
ORG_STRUC.E_FinImage = 10E+9;
ORG_STRUC.fp_FinImage=[];

% ****************************** %
% ****************************** %
% *     fingerprints           * %
ORG_STRUC.RmaxFing = 10;
ORG_STRUC.deltaFing = 0.04;
ORG_STRUC.sigmaFing = 0.03;
ORG_STRUC.toleranceFing = 0.005;

% Erf table for fast fingerprint calculations
ORG_STRUC.erf_table = zeros(803,1);
for i = 1 : 803
    ORG_STRUC.erf_table(i) = erf((i-402)/100);
end

% for compatibility
ORG_STRUC.varcomp = 0;
ORG_STRUC.molecule= 0;
ORG_STRUC.submitCount = 0;
