function createORG_AbinitCode(inputFile)

% USPEX Version 9.3.0
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%[nothing, abinitioCode] = unix(['./getFromTo abinitioCode ENDabinit ' inputFile]);
%ORG_STRUC.abinitioCode = str2num(abinitioCode);
ORG_STRUC.abinitioCode = python_uspex(getPy, ['-f ' inputFile ' -b abinitioCode -e ENDabinit'], 1);

ORG_STRUC.restartFile  = createRestartFileString(ORG_STRUC.abinitioCode);

%[nothing, whichCluster] = unix(['./getStuff ' inputFile ' whichCluster 1']);
whichCluster = python_uspex(getPy, ['-f ' inputFile ' -b whichCluster -c 1']);
if ~isempty(whichCluster)
    if ~isempty(str2num(whichCluster)) % number
        ORG_STRUC.platform=abs(str2num(whichCluster));
    else
        whichCluster(end) = [];
        if strcmp(whichCluster, 'nonParallel')
            ORG_STRUC.platform = 0;
        elseif strcmp(whichCluster, 'CFN')
            ORG_STRUC.platform = 3;
        elseif strcmp(whichCluster, 'QSH')
            ORG_STRUC.platform = 4;
        elseif strcmp(whichCluster, 'QSH2')
            ORG_STRUC.platform = 5;
        elseif strcmp(whichCluster, 'xservDE')
            ORG_STRUC.platform = 6;
        elseif strcmp(whichCluster, 'MIPT')
            ORG_STRUC.platform = 7;
        end
    end
end
if ORG_STRUC.platform > 7
    disp(['you did not chose a valid platform. Program STOPS...']);
    fprintf(fp,'you did not chose a valid platform. Program STOPS...\n');
    quit;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%[nothing, Kresol] =unix (['./getFromTo KresolStart Kresolend ' inputFile]);
Kresol = python_uspex(getPy, ['-f ' inputFile ' -b KresolStart -e Kresolend']);
if isempty(Kresol) % default
    Ltmp = length(ORG_STRUC.abinitioCode);
    ORG_STRUC.Kresol = 0.2*ones(1, Ltmp);
    if Ltmp > 1
        for i = 1 : Ltmp
            ORG_STRUC.Kresol(i) = 0.2 - (i-1)*0.12/(Ltmp-1);
        end
    end
else
    ORG_STRUC.Kresol = str2num(Kresol);
end


%[nothing, remoteFolder1] = unix (['./getStuff ' inputFile ' remoteFolder 1']);
remoteFolder1 = python_uspex(getPy, ['-f ' inputFile ' -b remoteFolder -c 1']);
ORG_STRUC.remoteFolder = remoteFolder1(1:length(remoteFolder1)-1);

% number of processors involved in each calculation
%[nothing, numProcessors] = unix(['./getFromTo numProcessors EndProcessors ' inputFile]);
numProcessors = python_uspex(getPy, ['-f ' inputFile ' -b numProcessors -e EndProcessors']);
if isempty(numProcessors)
    ORG_STRUC.numProcessors = ones(1, length(ORG_STRUC.abinitioCode))*8; % default
else
    ORG_STRUC.numProcessors = str2num(numProcessors);
end

%[nothing, maxErrors] =unix (['./getStuff ' inputFile ' maxErrors 1']);
maxErrors = python_uspex(getPy, ['-f ' inputFile ' -b maxErrors -c 1']);
if ~isempty(maxErrors)
    ORG_STRUC.maxErrors=str2num(maxErrors);
end

%[nothing, numParallelCalcs] = unix (['./getStuff ' inputFile ' numParallelCalcs 1']);
numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if ~isempty(numParallelCalcs)
    if str2num(numParallelCalcs) >= 2
        ORG_STRUC.numParallelCalcs = 2;
    else
        ORG_STRUC.numParallelCalcs = 1  ;
    end
end

%%%%%%%%%%%%% Command stuff
if ORG_STRUC.platform==0
    %[nothing, commandExecutable] = unix (['./getFromTo commandExecutable EndExecutable ' inputFile]);
    commandExecutable = python_uspex(getPy, ['-f ' inputFile ' -b commandExecutable -e EndExecutable']);
    helper = double(commandExecutable);
    tempIND = find (helper==10);
    helpIND = zeros(length(tempIND)+1,1);
    helpIND(2:end)=tempIND';
    for listLoop = 1:length(helpIND)-1
        listCommand{listLoop}=commandExecutable(helpIND(listLoop)+1:helpIND(listLoop+1)-1);
    end
    
    for commandLoop = 1: length(ORG_STRUC.abinitioCode)
        if ORG_STRUC.abinitioCode(commandLoop) == 0
            tempCom{commandLoop} = '0';     % means we do not do any optimization at all! (used in order optimization)
        elseif listLoop < length(ORG_STRUC.abinitioCode)
            tempCom{commandLoop} = listCommand{1};
        else
            tempCom{commandLoop} = listCommand{commandLoop};
        end
    end
    ORG_STRUC.commandExecutable=tempCom;
end
%%%%%%%%%%%%% END command stuff



function restartFile  = createRestartFileString(abinitioCode)

if isempty( abinitioCode )
    disp('Error : USPEX dont know which type of abinit should run, please check your input file... ')
else
    switch abinitioCode(1)
        case 4
            restartFile = 'lammps.restart';
        case 7
            restartFile = 'cp2k.restart';
    end
end
            
    
