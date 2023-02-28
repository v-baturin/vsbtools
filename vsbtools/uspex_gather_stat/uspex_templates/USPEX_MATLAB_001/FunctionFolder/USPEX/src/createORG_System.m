function createORG_System(inputFile)

% USPEX Version 9.3.0
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
global ORG_STRUC

%%%%%%% GET the code and cluster.
getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Whether we collect the forces and atomic position during the structure relaxation
collectForces = python_uspex(getPy, ['-f ' inputFile ' -b collectForces -c 1']);
if ~isempty(collectForces)
    ORG_STRUC.collectForces = str2num(collectForces);
end


%%%% Whether we fix the rand generator seeds
fixRndSeed = python_uspex(getPy, ['-f ' inputFile ' -b fixRandSeed -c 1']);
if ~isempty(fixRndSeed)
    ORG_STRUC.fixRndSeed = str2num(fixRndSeed);
end

%%%%%%% optimization parameter; optimize by: enthalpy, volume, hardness, etc
%[nothing, optType] = unix(['./getStuff ' inputFile ' optType 1']);
optType = python_uspex(getPy, ['-f ' inputFile ' -b optType -c 1']);
if ~isempty(optType)
    if optType(1) == '-'  % non documented: '-' at start means 'reverse' optimization (minimization <=> maximization)
        ORG_STRUC.opt_sign = -1;
        optType(1) = [];
    end
    optType(end) = [];

    if ~isempty(str2num(optType)) % number
        ORG_STRUC.optType=abs(str2num(optType));
    else
        if strcmpi(optType, 'enthalpy')
            ORG_STRUC.optType = 1;
        elseif strcmpi(optType, 'volume')
            ORG_STRUC.optType = 2;
        elseif strcmpi(optType, 'hardness')
            ORG_STRUC.optType = 3;
        elseif strcmpi(optType, 'struc_order')
            ORG_STRUC.optType = 4;
        elseif strcmpi(optType, 'density')
            ORG_STRUC.optType = 5;
        elseif strcmpi(optType, 'diel_sus')
            ORG_STRUC.optType = 6;
        elseif strcmpi(optType, 'gap')
            ORG_STRUC.optType = 7;
        elseif strcmpi(optType, 'diel_gap')
            ORG_STRUC.optType = 8;
        elseif strcmpi(optType, 'mag_moment')
            ORG_STRUC.optType = 9;
        elseif strcmpi(optType, 'quasientropy')
            ORG_STRUC.optType = 10;
	elseif strcmpi(optType, 'birefringence')
	    ORG_STRUC.optType = 11;
        elseif strcmpi(optType, 'TE_property')
            ORG_STRUC.optType = 14;
            ORG_STRUC.TEparam = struct({'Tmax', 800.0, 'Tdelta', 50.0, 'efcut', 0.15, 'goal', 'ZT'});
	elseif strcmpi(optType, 'Fphon')
	    ORG_STRUC.optType = 17;
	    ORG_STRUC.SCPHParam = struct('DISP',{},'Temp',{},'SCPHloop',{},'MPgrid',{},'Npoints',{},'DOSinputs',{});
        end
    end

    optType_OK = [1:11, 14,17, 1101:1111]; %% registered TE_property,Fphon
    if isempty( find(ORG_STRUC.optType == optType_OK) )
        fprintf(fp, 'you did not chose a valid optType. Program STOPS...\n');
        disp(['you did not chose a valid optType. Program STOPS...']);
        quit;
    end
end
% 1:'enthalpy'
% 2:'volume'
% 3:'hardness'
% 4:'struc_order'
% 5:'aver_dist'
% 6:'diel_sus'
% 7:'gap'
% 8:'diel_gap'
% 9:'mag_moment'
%10:'quasientropy'
%14:'TE_property'
%17:Free_energy 

%%%%%CalcType
%   s/S+3 digits (dimension: 0-3; molecule: 0/1; varcomp: 0/1)
%[nothing, calculationType] = unix(['./getStuff ' inputFile ' calculationType 1']);
calculationType = python_uspex(getPy, ['-f ' inputFile ' -b calculationType -c 1']);
if ~isempty(calculationType)
    %----------------- check Spin option -----------------%
    isSpinUp  = findstr(calculationType, 'S'); 
    isSpinLow = findstr(calculationType, 's');
    if ~isempty(isSpinUp) || ~isempty(isSpinLow)
        ORG_STRUC.spin = 1;
        ORG_STRUC.fracSpin = 0.1;
    end
    calculationType(isSpinUp)  = [];
    calculationType(isSpinLow) = [];
    %------------------------------------------------------%
    if calculationType(1) == '-'  % 2D
        calculationType(1) = [];
        ORG_STRUC.dimension = -1*str2num(calculationType(1));
    else
        ORG_STRUC.dimension = str2num(calculationType(1));
    end

    ORG_STRUC.molecule  = str2num(calculationType(2));
    ORG_STRUC.varcomp   = str2num(calculationType(3));
end

%%%%%Plugin Code
pluginType = python_uspex(getPy, ['-f ' inputFile ' -b pluginType -c 1']);
if ~isempty(pluginType)
    ORG_STRUC.pluginType = str2num(pluginType);
    if 1==ORG_STRUC.pluginType
        ORG_STRUC.startNextGen = 0;
    end
else
    ORG_STRUC.pluginType = 0;
end

% number of ions of each type
%[nothing, numIons] = unix(['./getFromTo numSpeci EndNumSpeci ' inputFile]);
%ORG_STRUC.numIons = str2num(numIons);
numIons = python_uspex(getPy, ['-f ' inputFile ' -b numSpeci -e EndNumSpeci'], 1);
ORG_STRUC.numIons = numIons;
if (ORG_STRUC.dimension==0) && (ORG_STRUC.varcomp==1) %001 mode
    ORG_STRUC.numIons = numIons';
end
if ORG_STRUC.molecule == 1
    ORG_STRUC.numMols = ORG_STRUC.numIons;
end
% types of ions; can be number, short name of full name. USPEX will use the numbers only
%[nothing, atomType] = unix(['./getFromTo atomType EndAtomType ' inputFile]);
atomType = python_uspex(getPy, ['-f ' inputFile ' -b atomType -e EndAtomType']);
N_type   = size(ORG_STRUC.numIons, 2);
atomType1 = GetElement(N_type, atomType);

for i = 1 : length(atomType1)
    if atomType1(i) == 0
        atomType1(i) = [];
    end
end
ORG_STRUC.atomType = atomType1;

% valences for each type of atoms
%[nothing, valences] = unix(['./getFromTo valences endValences ' inputFile]);
%ORG_STRUC.valences = str2num(valences);
valences = python_uspex(getPy, ['-f ' inputFile ' -b valences -e EndValences']);
if isempty(valences) % default
    ORG_STRUC.valences = zeros(1, length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        ORG_STRUC.valences(i) = str2num(valence(ORG_STRUC.atomType(i)));
    end
else
    ORG_STRUC.valences = str2num(valences);
end

% number of valence electrons for each type of atoms
%[nothing, NvalElectrons] = unix(['./getFromTo valenceElectr endValenceElectr ' inputFile]);
NvalElectrons = python_uspex(getPy, ['-f ' inputFile ' -b valenceElectr -e EndValenceElectr']);
if isempty(NvalElectrons) % default
    ORG_STRUC.NvalElectrons = zeros(1, length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        ORG_STRUC.NvalElectrons(i) = str2num(valenceElectronsNumber(ORG_STRUC.atomType(i)));
    end
else
    ORG_STRUC.NvalElectrons = str2num(NvalElectrons);
end

minAt = python_uspex(getPy, ['-f ' inputFile ' -b minAt -c 1']);
maxAt = python_uspex(getPy, ['-f ' inputFile ' -b maxAt -c 1']);
if ORG_STRUC.varcomp == 0 && (abs(ORG_STRUC.dimension) == 3 || ORG_STRUC.dimension == -2)
    if isempty(minAt) && isempty(maxAt)
        % Do nothing
    elseif   ~isempty(minAt) && ~isempty(maxAt)
        ORG_STRUC.minAt = str2num(minAt);
        ORG_STRUC.maxAt = str2num(maxAt);
    elseif   isempty(minAt) || isempty(maxAt)
        status = 'Please specify BOTH maximum and minimum amount of atoms in the structure for variable composition calculation'
        quit;
    end
elseif ORG_STRUC.varcomp == 1 && ORG_STRUC.dimension == 3
    % minimum amount of atoms/cell for first generation for varcomp
    %[nothing, minAt] = unix (['./getStuff ' inputFile ' minAt 1']);
    %minAt = python_uspex(getPy, ['-f ' inputFile ' -b minAt -c 1']);
    if isempty(minAt)
        status = 'Please specify the minimum amount of atoms in the structure for variable composition calculation'
        quit;
    end
    ORG_STRUC.minAt = str2num(minAt);
    % maximum amount of atoms/cell for first generation for varcomp
    %[nothing, maxAt] = unix (['./getStuff ' inputFile ' maxAt 1']);
    %maxAt = python_uspex(getPy, ['-f ' inputFile ' -b maxAt -c 1']);
    if isempty(maxAt)
        status = 'Please specify the maximum amount of atoms in the structure for variable composition calculation'
        quit;
    end
    ORG_STRUC.maxAt = str2num(maxAt);
end

%[nothing, ExternalPressure] = unix(['./getStuff ' inputFile ' ExternalPressure 1']);
%ORG_STRUC.ExternalPressure = str2num(ExternalPressure);
ExternalPressure = python_uspex(getPy, ['-f ' inputFile ' -b ExternalPressure -c 1']);
if ~isempty(ExternalPressure)
    ORG_STRUC.ExternalPressure = str2num(ExternalPressure);
else
    ORG_STRUC.ExternalPressure = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constraint
%[nothing, minVectorLength] = unix(['./getStuff ' inputFile ' minVectorLength 1']);
minVectorLength = python_uspex(getPy, ['-f ' inputFile ' -b minVectorLength -c 1']);
if isempty(minVectorLength)
    minVectorLength = 0;
    Vector = zeros(1, length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        Vector(i) = 1.8 * str2num(covalentRadius(ceil(ORG_STRUC.atomType(i))));
    end
    if ORG_STRUC.varcomp == 0
        minVectorLength = max(Vector);
    else
        minVectorLength = min(Vector);
    end
    ORG_STRUC.minVectorLength = minVectorLength; % default
else
    ORG_STRUC.minVectorLength = str2num(minVectorLength);
end

%[nothing, hardCore] = unix(['./getFromTo IonDistances EndDistances ' inputFile]);
hardCore = python_uspex(getPy, ['-f ' inputFile ' -b IonDistances -e EndDistances']);
if isempty(hardCore)
    tempFolder = [ORG_STRUC.homePath '/CalcFoldTemp'];
    codeFolder = [ORG_STRUC.USPEXPath '/FunctionFolder/Tool'];
    V_atom     = zeros(1, length(ORG_STRUC.atomType));
    hardCore   = zeros(length(ORG_STRUC.atomType), length(ORG_STRUC.atomType));
    for i = 1 : length(ORG_STRUC.atomType)
        V_atom(i) = calcDefaultVolume(1, ORG_STRUC.atomType(i), ORG_STRUC.ExternalPressure, 0, tempFolder, codeFolder);
    end
    for i = 1 : length(ORG_STRUC.atomType)
        for j = i:length(ORG_STRUC.atomType)
            if ORG_STRUC.molecule == 0
               hardCore(i,j) = min(0.22 * (V_atom(i) ^ (1/3) + V_atom(j) ^ (1/3)), 1.2);
            else
               hardCore(i,j) = 0.45 * (V_atom(i) ^ (1/3) + V_atom(j) ^ (1/3));
            end
        end
    end
else
    hardCore = str2num(hardCore);
end


%%%% check whether hard core radii or distance matrix were given and create matrix if necessary %%%%%%%
if size(hardCore, 1) == size(hardCore, 2)   %N*N matrix, read from INPUT.txt
    for i = 1 : length(ORG_STRUC.atomType)
        for j = i : length(ORG_STRUC.atomType)
            ORG_STRUC.minDistMatrice(i,j) = hardCore(i,j);
            ORG_STRUC.minDistMatrice(j,i) = hardCore(i,j);
        end
    end
%else %create matrix by default value
%    ORG_STRUC.minDistMatrice = zeros(length(hardCore));
%    for hardInd_1 = 1 : length(hardCore)
%        for hardInd_2 = hardInd_1 : length(hardCore)
%            ORG_STRUC.minDistMatrice(hardInd_1,hardInd_2) = hardCore(hardInd_1) + hardCore(hardInd_2);
%            ORG_STRUC.minDistMatrice(hardInd_2,hardInd_1) = hardCore(hardInd_1) + hardCore(hardInd_2);
%        end
%    end
end
%%%%%%%%%%%%%% END Constraint %%%%%%%%%%%%%%%%%%

%%%%%% Good Bonds  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bonds with nu higher than this value are always included into hardness formula
ORG_STRUC.goodBonds = zeros(length(ORG_STRUC.atomType));
%[nothing, goodBonds] = unix(['./getFromTo goodBonds EndGoodBonds ' inputFile]);
goodBonds = python_uspex(getPy, ['-f ' inputFile ' -b goodBonds -e EndGoodBonds']);
if isempty(goodBonds)
    gB = calcDefaultGoodBonds(ORG_STRUC.atomType); % default
else
    gB = str2num(goodBonds);
end
if size(gB,1) == 1        % just 1 value present or 1 type of atoms
    for i = 1 : length(ORG_STRUC.atomType)
        for j = 1 : length(ORG_STRUC.atomType)
            ORG_STRUC.goodBonds(i,j) = gB;
        end
    end
else
    for i = 1 : length(ORG_STRUC.atomType)
        for j = i : length(ORG_STRUC.atomType)
            ORG_STRUC.goodBonds(i,j) = gB(i,j);
            ORG_STRUC.goodBonds(j,i) = gB(i,j);
        end
    end
end

%%%%%%%%% END GoodBonds %%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%% magnetic type ratio

% magnetic type ratio
%[nothing, magRatioVal] = unix(['./getFromTo magRatio EndMagRatio ' inputFile]);
magRatioVal = python_uspex(getPy, ['-f ' inputFile ' -b magRatio -e EndMagRatio']);
if isempty(magRatioVal)
    ORG_STRUC.magRatio(1:7) = [0.1, 0.9/4, 0.9/4, 0.9/4, 0.9/4, 0, 0];
else
    ORG_STRUC.magRatio = str2num(magRatioVal);
    ORG_STRUC.magRatio = abs(ORG_STRUC.magRatio(1:7))/sum(abs(ORG_STRUC.magRatio));
end

% ldaU 
ldaU = python_uspex(getPy, ['-f ' inputFile ' -b ldaU -e EndLdaU']);
if isempty(ldaU)
    ORG_STRUC.ldaU(1 : length(ORG_STRUC.atomType)) = 0;
else
    ORG_STRUC.ldaU = str2num(ldaU);
    if size(ORG_STRUC.ldaU, 1) == 1;
        ORG_STRUC.ldaU(2, :) = 0;
    end
end

%% Thermoelectric property parameters
if 14 == ORG_STRUC.optType
  TE_goal = python_uspex(getPy, ['-f ' inputFile ' -b TE_goal -c 1']);
  TE_goal = strtrim(TE_goal);
  if ~isempty(TE_goal)
    ORG_STRUC.TEparam.goal = TE_goal;
  end
  input_fields = {'TE_threshold', 'TE_T_interest', 'BoltzTraP_T_max', 'BoltzTraP_T_delta', 'BoltzTraP_efcut'};
  struct_fields = {'threshold', 'Tinterest', 'Tmax', 'Tdelta', 'efcut'};
  for ii = 1:length(input_fields)
    field_name = input_fields{ii};
    value = python_uspex(getPy, ['-f ' inputFile ' -b ' field_name ' -c 1']);
    value = strtrim(value);
    if ~isempty(value)
      ORG_STRUC.TEparam.(struct_fields{ii}) = str2num(value);
    end
  end
  % ORG_STRUC.TEparam %% For debug
end
%SCAILD 
if 17 == ORG_STRUC.optType
       	DISP = python_uspex(getPy, ['-f ' inputFile ' -b DISP -c 1']);
	if isempty(DISP)
	       	ORG_STRUC.SCPHParam.DISP = 40;
	else
		ORG_STRUC.SCPHParam.DISP = str2num(DISP);
	 end
	 Temp = python_uspex(getPy, ['-f ' inputFile ' -b Temperature -c 1']);
	if isempty(Temp)
		status ='Please specify Temperature for Finite Temperature prediction';
	       	quit;
	 else
	       	ORG_STRUC.SCPHParam.Temp = str2num(Temp);
	end
	SCPHloop = python_uspex(getPy, ['-f ' inputFile ' -b SCPHloop -c 1']);
	if isempty(SCPHloop)
	       	ORG_STRUC.SCPHParam.SCPHloop = 30;%default
        else
		ORG_STRUC.SCPHParam.SCPHloop = str2num(SCPHloop);
	 end
	 Npoints = python_uspex(getPy, ['-f ' inputFile ' -b Npoints -c 1']);
	 if isempty(Npoints)
		 ORG_STRUC.SCPHParam.Npoints = 60;
	else
		ORG_STRUC.SCPHParam.Npoints = str2num(Npoints);
	end
	MPgrid = python_uspex(getPy, ['-f ' inputFile ' -b MPgrid -c 1']); 
	if isempty(MPgrid) 
		ORG_STRUC.SCPHParam.MPgrid = 60;
	else 
		ORG_STRUC.SCPHParam.MPgrid = str2num(MPgrid);
	 end 
	 Dosinputs = python_uspex(getPy, ['-f' inputFile ' -b Dosinputs -e EndDosinputs']);
	 ORG_STRUC.SCPHParam.Dosinputs = str2num(Dosinputs);
end

%% Machine learning based population prefiltering
mlPeerFilter = python_uspex(getPy, ['-f ' inputFile ' -b mlPeerFilter -c 1']);
if isempty(mlPeerFilter)
    ORG_STRUC.mlPeerFilter = 0;
else
    ORG_STRUC.mlPeerFilter = str2num(mlPeerFilter);
end
mlMinPopulation = python_uspex(getPy, ['-f ' inputFile ' -b mlMinPopulation -c 1']);
if isempty(mlMinPopulation)
    ORG_STRUC.mlMinPopulation = 128;
else
    ORG_STRUC.mlMinPopulation = str2num(mlMinPopulation);
end
mlFilterRatio = python_uspex(getPy, ['-f ' inputFile ' -b mlFilterRatio -c 1']);
if isempty(mlFilterRatio)
    ORG_STRUC.mlFilterRatio = 5;
else
    ORG_STRUC.mlFilterRatio = str2num(mlFilterRatio);
end

%%%%%%%%%%% Restart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if pickUpYN ~= 0 , then a previous calculation will be continued
%[nothing, pickUpYN] = unix (['./getStuff ' inputFile ' pickUpYN 1']);
%pickUpYN = python_uspex(getPy, ['-f ' inputFile ' -b pickUpYN -c 1']);
%if ~isempty(pickUpYN)
%    ORG_STRUC.pickUpYN = str2num(pickUpYN);
%end

% from which generation shall the last calculation be continued? If = 0, then the last available generation is taken
%[nothing, pickUpGen] = unix (['./getStuff ' inputFile ' pickUpGen 1']);
pickUpGen = python_uspex(getPy, ['-f ' inputFile ' -b pickUpGen -c 1']);
if ~isempty(pickUpGen)
    ORG_STRUC.pickUpGen = str2num(pickUpGen);
    if ORG_STRUC.pickUpGen > 0
        ORG_STRUC.pickUpYN  = 1;
    end
end
if ORG_STRUC.pickUpYN == 0
    ORG_STRUC.pickUpGen = 1; % needed for BestEnthalpies plot
end
% number of the results folder to be used. If = 0 , then the highest
% existing number is taken
%[nothing, pickUpFolder] =unix (['./getStuff ' inputFile ' pickUpFolder 1']);
pickUpFolder = python_uspex(getPy, ['-f ' inputFile ' -b pickUpFolder -c 1']);
if ~isempty(pickUpFolder)
    ORG_STRUC.pickUpFolder = str2num(pickUpFolder);
end
%%%%%%%%%%End Restart %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
