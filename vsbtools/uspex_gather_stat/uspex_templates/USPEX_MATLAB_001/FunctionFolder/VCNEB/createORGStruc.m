function createORGStruc(inputFile)


global ORG_STRUC

%--------------------------------------------------------------------------
%   Load Default ORG_STRUC for VCNEB

memStruct = NEB_creatORGDefault();

ORG_STRUC.inputFile = inputFile;

%--------------------------------------------------------------------------


% %%%%%%%% RESULTS
existfold = 1;
folderNum = 0;
while ~isempty(existfold)
    folderNum = folderNum + 1;
    ORG_STRUC.resFolder = ['results' num2str(folderNum)];
    [nothing1,existfold] = mkdir(ORG_STRUC.resFolder);
end
mkdir( [ORG_STRUC.resFolder '/PATH'] );
mkdir( [ORG_STRUC.resFolder '/STEP'] );
mkdir( [ORG_STRUC.resFolder '/AuxiliaryFiles'] );


fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'w');


% Create a header automatically:
description    = 'Variable-Cell Nudged-Elastic-Band Method for Phase Transition Investigation';
funcFold       = [ORG_STRUC.USPEXPath '/FunctionFolder'];
formatted_rows = createHeader('USPEX_pic.txt', description, funcFold);

for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end


% Cite:
text = {'Please cite the following suggested papers',       ...
        'when you publish the results obtained from USPEX:' ...
        };
formatted_rows = createHeader_wrap(text);
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end

text = {'Qian G.R., Dong X., Zhou X.F., Tian Y.J., Oganov A.R., Wang H.T. (2013)', ...
        'Variable cell nudged elastic band method for studying solid-solid',  ...
        'structural phase transitions.', ...
        'Comp. Phys. Comm., 184, 2111-2118', ...
        '', ...
        'X. Dong, X.F. Zhou, G.R. Qian, Z.S. Zhao, Y.J. Tian, H.T. Wang (2013)', ...
        'An ab initio study on the transition paths from graphite', ...
        'to diamond under pressure', ...
        'J. Phys.: Condens. Matter 25, 145402' ...
        };
formatted_rows = createHeader_wrap(text, 'left');
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end

text = {datestr(now)};
formatted_rows = createHeader_wrap(text);
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end
fprintf(fp,'\n');

text = {'', ...
        'Please see VCNEBReports for brief results', ...
        '', ...
        };
formatted_rows = createHeader_wrap(text, 'center');
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end
fprintf(fp,'\n');
fclose(fp);


%-----------------------------------------------------------------------------------------------------
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%%%%% vcnebType
%3 digits  (VCNEB/Structural Relaxation: 0/1; Variable-Image: 0/1;  Variable-Spring Constant: 0/1   )

%[nothing, vcnebType] = unix(['./getStuff ' inputFile ' vcnebType 1']);
vcnebType = python_uspex(getPy, ['-f ' inputFile ' -b vcnebType -c 1']);
if ~isempty(vcnebType)
    ORG_STRUC.CalcType = str2num(vcnebType(1));
    ORG_STRUC.optVarImage  = str2num(vcnebType(2));
    ORG_STRUC.optVarK  = str2num(vcnebType(3));
end

for i = 1:length( memStruct.numeral )
    %    [nothing, getResult] = unix(['./getStuff ' inputFile ' ' memStruct.numeral{i} ' 1']);
    getResult = python_uspex(getPy, ['-f ' inputFile ' -b '  memStruct.numeral{i} ' -c 1']);
    if ~isempty(getResult)
        ORG_STRUC = setfield(ORG_STRUC, memStruct.numeral{i}, str2num(getResult) );
    end
end


for i = 1:length( memStruct.string )
    %    [nothing, getResult] = unix(['./getStuff ' inputFile ' ' memStruct.string{i} ' 1']);
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
% 
% ExternalPressure = python_uspex(getPy, ['-f ' inputFile ' -b ExternalPressure -e EndPressure']);
% ORG_STRUC.ExternalPressure = str2num(ExternalPressure);%/10;
%-----------------------------------------------------------------------------

% Read Abinit Code parameters
createORG_AbinitCode(inputFile);

% Read Fingerprint parameters
createORG_Fingerprint(inputFile);

%--------------------------------------------------------------
% types of ions; can be number, short name of full name. USPEX will use the numbers only
%[nothing, atomType] = unix(['./getFromTo atomType EndAtomType ' inputFile]);
atomType = python_uspex(getPy, ['-f ' inputFile ' -b atomType -e EndAtomType']);
atomType1 = zeros(1,size(ORG_STRUC.numIons,2));
c1 = findstr(atomType, ' ');
c = sort(str2num(['0 ' num2str(c1)]));
c(end+1) = length(atomType) + 1;
ind1 = 1;
for i = 2 : length(c)
    if c(i-1)+1 > c(i)-1
        continue
    end
    tmp = atomType(c(i-1)+1 : c(i)-1);
    tmp = strtrim(tmp);
    if ~isempty(str2num(tmp)) % number
        atomType1(ind1) = str2num(tmp);
    else
        if strcmp(tmp,'H.5')
            atomType1(ind1)=0.5;
        elseif strcmp(tmp,'H.75')
            atomType1(ind1)=0.75;
        elseif strcmp(tmp,'H1.25')
            atomType1(ind1)=1.25;
        elseif strcmp(tmp,'H1.5')
            atomType1(ind1)=1.5;
        else
            for j = 1 : 105
                if strcmp(lower(tmp), lower(elementFullName(j))) || strcmp(lower(tmp), lower(megaDoof(j)))
                    atomType1(ind1) = j;
                    break;
                end
            end
        end
    end
    ind1 = ind1 + 1;
end


for i=1:length(atomType1)
    if atomType1(i)==0
        atomType1(i)=[];
    end
end
ORG_STRUC.atomType = atomType1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ORG_STRUC.atomSymbol='';
for i = 1:length(ORG_STRUC.numIons)
    ORG_STRUC.atomSymbol{i}=  megaDoof( ORG_STRUC.atomType(i) );
end
%disp('Atomic Info : ')
% ORG_STRUC.atomType
%for i = 1:length(ORG_STRUC.numIons)
%    ORG_STRUC.atomSymbol{i}
%end
% ionChange is just a handy little info thing
for ind = 1:length(ORG_STRUC.numIons)
    ORG_STRUC.ionChange(ind) = sum(ORG_STRUC.numIons(1:ind));
end

ORG_STRUC.ionCh = zeros(1,length(ORG_STRUC.ionChange)+1);
for i = 2 : length(ORG_STRUC.ionChange)+1
    ORG_STRUC.ionCh(i) = ORG_STRUC.ionChange(i-1);
end
%-------------------------------------------------------------


% Space group determination tolerance
%[nothing, SGtolerance] = unix(['./getStuff ' inputFile ' SymTolerance 1']);
SGtolerance = python_uspex(getPy, ['-f ' inputFile ' -b SymTolerance -c 1']);
if ~isempty(SGtolerance)
    if strcmp(lower(SGtolerance), 'high')
        ORG_STRUC.SGtolerance = 0.05;
    elseif strcmp(lower(SGtolerance), 'medium')
        ORG_STRUC.SGtolerance = 0.10;
    elseif strcmp(lower(SGtolerance), 'low')
        ORG_STRUC.SGtolerance = 0.20;
    else % number!
        ORG_STRUC.SGtolerance = str2num(SGtolerance);
    end
end
%ORG_STRUC


% Defalut optimize algorithm in VCNEB
if  ORG_STRUC.optimizerType==0
    if  ORG_STRUC.CalcType==1      % SD for VCNEB
        ORG_STRUC.optimizerType=1;
    elseif ORG_STRUC.CalcType==2   % FIRE for relaxation
        ORG_STRUC.optimizerType=2;
    end
end


% the K-max,K-min and K-constant
if ORG_STRUC.optVarK==0
    ORG_STRUC.K_min=ORG_STRUC.Kconstant;
    ORG_STRUC.K_max=ORG_STRUC.Kconstant;
end


%ORG_STRUC

%---------------------------------------------------------------
NEB_constraint();
