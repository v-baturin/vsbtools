function createORG_AbinitCode(inputFile)

global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

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


maxErrors = python_uspex(getPy, ['-f ' inputFile ' -b maxErrors -c 1']);
if ~isempty(maxErrors)
   ORG_STRUC.maxErrors = str2num(maxErrors);
end

% ORG_STRC.conv_till = how many steps are optimized at constant volume; rest - 'final' optimization (not used in ranking)
% format - 1 1 1 1 (1 1) means that 2 last steps are 'final' optimization steps
abinitioCode = python_uspex(getPy, ['-f ' inputFile ' -b abinitioCode -e ENDabinit']);
if isempty(findstr(abinitioCode, '('))
  ORG_STRUC.abinitioCode = str2num(abinitioCode);
  ORG_STRUC.conv_till = length(ORG_STRUC.abinitioCode);
else
  ORG_STRUC.conv_till = length(str2num(abinitioCode(1:findstr(abinitioCode,'(')-1)));
  tmp = strrep(abinitioCode, '(', ' ');
  tmp1 = strrep(tmp, ')', ' ');
  ORG_STRUC.abinitioCode = str2num(tmp1);
end

% %%%%%%%% RESULTS
ExternalPressure = python_uspex(getPy, ['-f ' inputFile ' -b ExternalPressure -c 1']);
if isempty(ExternalPressure)
   ExternalPressure = '0'; % default
end
ORG_STRUC.ExternalPressure = str2num(ExternalPressure);


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

%%%%%%%%%%%%% Command stuff
if ORG_STRUC.platform==0
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
FullRelax = python_uspex(getPy, ['-f ' inputFile ' -b FullRelax -c 1']);
if isempty(FullRelax)
  FullRelax = '2'; % default
end
ORG_STRUC.FullRelax = str2num(FullRelax);

remote = python_uspex(getPy, ['-f ' inputFile ' -b remoteRegime -c 1']);
if isempty(remote)
  remote = '0'; % default
end
ORG_STRUC.remote = str2num(remote);

remoteFolder1 = python_uspex(getPy, ['-f ' inputFile ' -b remoteFolder -c 1']);
ORG_STRUC.remoteFolder = remoteFolder1(1:length(remoteFolder1)-1);

numParallelCalcs = python_uspex(getPy, ['-f ' inputFile ' -b numParallelCalcs -c 1']);
if isempty(numParallelCalcs)
  numParallelCalcs = '1'; % default
end
ORG_STRUC.numParallelCalcs = str2num(numParallelCalcs);

% number of processors involved in each calculation
numProcessors = python_uspex(getPy, ['-f ' inputFile ' -b numProcessors -e EndProcessors']);
if isempty(numProcessors)
  ORG_STRUC.numProcessors = ones(1, length(ORG_STRUC.abinitioCode))*8; % default
else
  ORG_STRUC.numProcessors = str2num(numProcessors);
end
