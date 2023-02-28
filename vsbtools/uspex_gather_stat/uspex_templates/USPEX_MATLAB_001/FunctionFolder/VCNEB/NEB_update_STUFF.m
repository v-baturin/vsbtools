function NEB_update_STUFF(inputFile)

global ORG_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Some parameters can be updated from INPUT.txt %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------------
memStruct.numeral = {
    'optVarImage';   'optVarK';
    'optimizerType'; 'optMethodCIDI'; 'optFreezing';
    'whetherConstraint';
    'dt'; 'K_min'; 'K_max'; 'Kconstant';
    'VarPathLength';
    'startCIDIStep';
    'ConvThreshold';
    %    'numSteps';
    'FormatType'; 'PrintStep';
    'SGtolerance';
    'numProcessors';     'numParallelCalcs';
    };

memStruct.array = {
    'whichCI';      'whichDI';
    'pickupImages';
    'Kresol';
    };

memStruct.string = {};

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];


vcnebType = python_uspex(getPy, ['-f ' inputFile ' -b vcnebType -c 1']);
if ~isempty(vcnebType)
    if str2num(vcnebType(1)) ~= ORG_STRUC.CalcType
        USPEXmessage(201,'',0);
    end
    %    ORG_STRUC.CalcType = str2num(vcnebType(1));
    ORG_STRUC.optVarImage  = str2num(vcnebType(2));
    ORG_STRUC.optVarK      = str2num(vcnebType(3));
end

%----------------------------------------------------------------------------
for i = 1:length( memStruct.numeral )
    getResult = python_uspex(getPy, ['-f ' inputFile ' -b '  memStruct.numeral{i} ' -c 1']);
    if ~isempty(getResult)
        ORG_STRUC = setfield(ORG_STRUC, memStruct.numeral{i}, str2num(getResult) );
    end
end

%----------------------------------------------------------------------------
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
    getResult = python_uspex(getPy, ['-f ' inputFile ' -b ' iniStr ' -e ' endStr ]);
    if ~isempty(getResult)
        ORG_STRUC = setfield(ORG_STRUC, memStruct.array{i}, str2num(getResult) );
    end
end


%----------------------------------------------------------------------------
if isempty(ORG_STRUC.numIons) && ~isempty(ORG_STRUC.numSpecies)
    ORG_STRUC.numIons = ORG_STRUC.numSpecies;  % compatible with USPEX, no numIons in INPUT.txt
end

% Defalut optimize algorithm in VCNEB
if  ORG_STRUC.optimizerType==0
    if  ORG_STRUC.CalcType==1      % SD for VCNEB
        ORG_STRUC.optimizerType=1;
    elseif ORG_STRUC.CalcType==2   % FIRE for relaxation
        ORG_STRUC.optimizerType=2;
    end
end


% the K-max,K-min and K-constant
if ORG_STRUC.optVarK==1
    ORG_STRUC.K_min=ORG_STRUC.Kconstant;
    ORG_STRUC.K_max=ORG_STRUC.Kconstant;
end




%---------------------------------------------------------------
NEB_constraint();
CreateCalcFolder();


safesave ('Current_ORG.mat', ORG_STRUC)

