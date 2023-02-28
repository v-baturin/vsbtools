function MINHOP_CreateORGStruc(inputFile)

global ORG_STRUC

ORG_STRUC = struct('dimension',{},'varcomp',{},'molecule',{},...       $system
    'minDistMatrice',{}, 'numIons',{},...   %atoms
    'atomType', {}, 'lattice', {}, 'coordinates', {}, ...              %atoms
    'valences', {}, 'NvalElectrons', {}, 'goodBonds',{}, 'checkConnectivity', {},...
    'numGenerations', {}, 'populationSize', {}, ...                
    'GaussianWidth', {}, 'GaussianHeight', {},...                      %metadynamics 
    'specificFolder',{}, 'homePath',{}, 'resFolder',{}, ...  %path
    'pickUpYN', {}, 'pickUpGen', {}, 'pickUpFolder',{}, 'pickedUP',{},...            %pickup
    'RmaxFing', {}, 'deltaFing', {}, 'sigmaFing', {}, 'doFing', {}, ... %finggerprint
    'toleranceFing', {}, 'erf_table', {}, 'weight', {}, 'log_file', {},...         %I-O
    'commandExecutable', {}, 'Kresol', {}, 'abinitioCode', {}, ...      %ab-init calculation
    'platform', {}, 'remoteFolder',{}, 'SGtolerance',{}, ...
    'maxErrors', {}, 'numProcessors', {},...
    'numParallelCalcs',{}, 'doSpaceGroup',{},'ExternalPressure',{}, ... 
    'Softmode_Fre', {}, 'Softmode_num', {},'Softmode_eignFre', {},...   %softmode
    'Softmode_eignVec', {}, 'Softmode_eignSupercells',{}, ...           %softmode
    'FullRelax', {},  'coorMutationDegree',{},'conv_till', {}, ...       %others
    'maxAt',{}, 'correctionAngle',{},  ...  %basics
    'basicStructureNumber', {}, 'useBasicCell', {},...                  %basics
    'maxVectorLength', {}, 'repeatForStatistics', {},...                                           %basics
    'maxIncrease', {}, 'minVectorLength', {}, 'correctionLength', {},...
    'constLattice', {}, 'minAngle' , {}, 'minDiagAngle' ,{}, 'mutationRate', {});  %nnn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply new path format:

[ORG_STRUC(1).homePath, ORG_STRUC(1).USPEXPath] = workingPath(); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ORG_STRUC.dimension = 3;
ORG_STRUC.molecule  = 0;
ORG_STRUC.varcomp   = 1;
ORG_STRUC.specificFolder = 'Specific';
ORG_STRUC.log_file = 'OUTPUT.txt';
ORG_STRUC.basicStructureNumber = 0;
ORG_STRUC.checkConnectivity = 1; % 
ORG_STRUC.maxErrors = 2;
ORG_STRUC.useBasicCell = 5;
ORG_STRUC.correctionAngle = 0;
ORG_STRUC.correctionLength = 0;
ORG_STRUC.repeatForStatistics = 1;
% ****************************** %
% ****************************** %
% *     fingerprints           * %  
% ****************************** %
ORG_STRUC.erf_table = zeros(803,1);
for i = 1 : 803
  ORG_STRUC.erf_table(i) = erf((i-402)/100);
end
ORG_STRUC.RmaxFing = 10;
ORG_STRUC.deltaFing = 0.08;
ORG_STRUC.sigmaFing = 0.03;
ORG_STRUC.toleranceFing = 0.008;
ORG_STRUC.toleranceBestHM = 0.02;
ORG_STRUC.constLattice = 1;
ORG_STRUC.mutationRate    =  0.05;
ORG_STRUC.minAngle = 55; %%ezafe shod
ORG_STRUC.minDiagAngle = 30;

ORG_STRUC.Softmode_num = 0;
createORG_System(inputFile);
createORG_AbinitCode(inputFile);
createORG_Symmetry(inputFile);
createORG_Fingerprint(inputFile);
createORG_MINHOP(inputFile);

% %%%%%%%% RESULTS
existfold = 1;
folderNum = 0;
while ~isempty(existfold)
    folderNum = folderNum + 1;
    ORG_STRUC.resFolder = ['results' num2str(folderNum)];
    [nothing1,existfold] = mkdir(ORG_STRUC.resFolder);
end

%%%%%%%read initial structure%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeORG();
