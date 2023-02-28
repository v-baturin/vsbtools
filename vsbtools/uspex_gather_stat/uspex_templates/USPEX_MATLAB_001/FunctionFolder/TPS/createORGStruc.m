function createORGStruc(inputFile)


global ORG_STRUC

%---------------------------------------------------
%   Load Default ORG_STRUC for VCNEB

memStruct = TPS_creatORGDefault();

%---------------------------------------------------
ORG_STRUC.submitCount = 0;

%--------------------------------------------------------------------------

% %%%%%%%% RESULTS
existfold = 1;
folderNum = 0;
while ~isempty(existfold)
    folderNum = folderNum + 1;
    ORG_STRUC.resFolder = ['results' num2str(folderNum)];
    [nothing1,existfold] = mkdir(ORG_STRUC.resFolder);
end


fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'w');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp,'|                                                                                 |\n');
fprintf(fp,'|                       __  __ _____  ____   ______ _  __                         |\n');
fprintf(fp,'|                      / / / // ___/ / __ \\ / ____/| |/ /                         |\n');
fprintf(fp,'|                     / / / / \\__ \\ / /_/ // __/   |   /                          |\n');
fprintf(fp,'|                    / /_/ / ___/ // ____// /___  /   |                           |\n');
fprintf(fp,'|                    \\____/ /____//_/    /_____/ /_/|_|                           |\n');
fprintf(fp,'|                                                                                 |\n');
fprintf(fp,'|                       Version 10.0.0 (12/05/2013)                               |\n'); 
fprintf(fp,'|                                                                                 |\n');
fprintf(fp,'|                     Transition Path Sampling Method                             |\n');
fprintf(fp,'|                more info at  http://uspex.stonyrbook.edu/                       |\n');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp,'|  Please cite the following suggested literatures                                |\n');
fprintf(fp,'|  when you publish the results obtained from USPEX                               |\n');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp,'|  Dellago, C., Bolhuis, P. G., Csajka, F. S. & Chandler, D.  (1998)              |\n');
fprintf(fp,'|  "Transition path sampling and the calculation of rate constants."              |\n');
fprintf(fp,'|  J. Chem. Phys. 108, 1964-1978                                                  |\n');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp,'|  Bolhuis, P. G., Chandler, D., Dellago, C. & Geissler, P. L.(2002)              |\n');
fprintf(fp,'|  "Transition path sampling:                                                     |\n');
fprintf(fp,'|               Throwing ropes over rough mountain passes, in the dark."          |\n');
fprintf(fp,'|  Annu. Rev. Phys. Chem. 53, 291-318                                             |\n');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp,'|  S. E. Boulfelfel, A. R. Oganov & S. Leoni(2012)                                |\n');
fprintf(fp,'|  "Understanding the nature of "superhard graphite"."                            |\n');
fprintf(fp,'|  Scientific Reports 2, 471                                                      |\n');
fprintf(fp,'*---------------------------------------------------------------------------------*\n');
fprintf(fp, ' -                        %20s                                   -\n', date() );
fprintf(fp,'\n');

fprintf(fp,'Iteration Calc1( op, Suc)->Calc2( op, Suc) HowCome     Amplitude    Shifter:   dH   -  Temp\n' );
fprintf(fp,'                                                     (A2B)   (B2A)  (accept)(Kcal/mol) (K)\n' );

fclose(fp);


%-----------------------------------------------------------------------------------------------------
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];


for i = 1:length( memStruct.numeral )
    getResult = python_uspex(getPy, ['-f ' inputFile ' -b '  memStruct.numeral{i} ' -c 1']);
    if ~isempty(getResult)
        ORG_STRUC = setfield(ORG_STRUC, memStruct.numeral{i}, str2num(getResult) );
    end
end


for i = 1:length( memStruct.string )
    getResult = python_uspex(getPy, ['-f ' inputFile ' -b '  memStruct.string{i} ' -c 1']);
    if ~isempty(getResult)
        getResult(end)=[];
        ORG_STRUC = setfield(ORG_STRUC, memStruct.string{i}, getResult );
    end
end

%----------------------------------------------------------------------------
for i = 1:length(  memStruct.array )
    iniStr = memStruct.array{i};
    endStr = [ 'End' upper(iniStr(1)) iniStr(2:end)];
    %    [nothing, getResult] = unix(['./getFromTo ' iniStr ' ' endStr ' ' inputFile  ]);
    getResult = python_uspex(getPy, ['-f ' inputFile ' -b ' iniStr ' -e ' endStr ]);
    if ~isempty(getResult)
        ORG_STRUC = setfield(ORG_STRUC, memStruct.array{i}, str2num(getResult) );
    end
end


%----------------------------------------------------------------------------
if isempty(ORG_STRUC.numIons) && ~isempty(ORG_STRUC.numSpecies)
    ORG_STRUC.numIons = ORG_STRUC.numSpecies;  % compatible with USPEX, no numIons in INPUT.txt
end

ExternalPressure = python_uspex(getPy, ['-f ' inputFile ' -b ExternalPressure -e EndPressure']);
ORG_STRUC.ExternalPressure = str2num(ExternalPressure);
%-----------------------------------------------------------------------------

% Read Abinit Code parameters
createORG_AbinitCode(inputFile);

%--------------------------------------------------------------
%atomType = python_uspex(getPy, ['-f ' inputFile ' -b atomType -e EndAtomType']);
%ORG_STRUC.atomType, ORG_STRUC.atomSymbol] = convertAtomType(atomType, ORG_STRUC.numIons);

%--------
ORG_STRUC.cmdCalcOp = python_uspex(getPy, ['-f ' inputFile ' -b cmdOrderParameter   -e EndCmdOrderParameter']);

ORG_STRUC.cmdCalcHT = python_uspex(getPy, ['-f ' inputFile ' -b cmdEnthalpyTempture -e EndCmdEnthalpyTempture']);


speciesSymbol = python_uspex(getPy, ['-f ' inputFile ' -b speciesSymbol -e EndSpeciesSymbol']);
ORG_STRUC.speciesSymbol = createSpeciesSymbol(speciesSymbol, ORG_STRUC.numSpecies);

if isempty( ORG_STRUC.mass )
    ORG_STRUC.mass = createSpeciesMass(ORG_STRUC.speciesSymbol);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ionChange is just a handy little info thing
%for ind = 1:size(ORG_STRUC.numIons,1)
%    ORG_STRUC.ionChange(ind) = sum(ORG_STRUC.numIons(1:ind));
%end

%ORG_STRUC.ionCh = zeros(1,length(ORG_STRUC.ionChange)+1);
%for i = 2 : length(ORG_STRUC.ionChange)+1
%    ORG_STRUC.ionCh(i) = ORG_STRUC.ionChange(i-1);
%end

%for i = 1:size(ORG_STRUC.numIons,1)
%    ORG_STRUC.mass(i) = 0;
%end

%-------------------------------------------------------------

% %%%%%%%% Create Results folders
%mkdir( [ORG_STRUC.resFolder '/Store_A'] );
%mkdir( [ORG_STRUC.resFolder '/Store_B'] );
mkdir( [ORG_STRUC.resFolder '/Iterations'] );
mkdir( [ORG_STRUC.resFolder '/AuxiliaryFiles'] );

%
%
%
function speciesSymbol = createSpeciesSymbol(speciesSymbolString, numSpecies) 


speciesSymbol = {}; 
c1 = findstr(speciesSymbolString, ' ');
c = sort(str2num(['0 ' num2str(c1)]));
c(end+1) = length(speciesSymbolString) + 1;

ind1 = 1;
for i = 2 : length(c)
    if c(i-1)+1 > c(i)-1
        continue
    end
    tmp = speciesSymbolString(c(i-1)+1 : c(i)-1);
    speciesSymbol{ind1} =  strtrim(tmp); 
    ind1 = ind1 + 1;
end



% inputCheck()  % very important here!! default fpTolerance shoule be set

%==========================================================================
%
%==========================================================================
function speciesMass = createSpeciesMass(speciesSymbol)

speciesMass = zeros(1, length(speciesSymbol) );
if isempty( speciesSymbol )
    error('No speciesSymbol has been set ...');
    %for i = 1:size(numIons,1)
    %    speciesSymbol(i,1:end);
    %end
else
    for i = 1:length( speciesSymbol )
        speciesStr  = speciesSymbol{i};
        [atomType, numIons]  = convertAtomType( speciesStr );
        for j = 1:length(atomType)
            speciesMass(i) = elementMass(atomType(j))*numIons(j)+speciesMass(i);
        end
    end
end


function [atomType, numIons]  = convertAtomType( speciesString )


speciesString=[ speciesString,'Z' ];
upChPos  = find(isstrprop(speciesString,'upper')>0);

atomType = zeros(1, length(upChPos)-1);
numIons  = zeros(1, length(upChPos)-1);

ind1 = 1;
for i = 2:length(upChPos)
    atomStr = speciesString(upChPos(i-1):upChPos(i)-1);
    digPos = find(isstrprop(atomStr,'digit')>0);
    if isempty(digPos)
        digPos=length(atomStr)+1;
        atomStr(end+1)='1';
    end
    atomTypeStr  = atomStr(1:digPos(1)-1);
    atomType(i-1)= elementID( atomTypeStr );
    numIons(i-1) = str2num( atomStr(digPos(1):end) );
end

