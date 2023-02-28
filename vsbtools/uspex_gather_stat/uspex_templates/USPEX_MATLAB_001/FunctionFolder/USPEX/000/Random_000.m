function Random_000(Ind_No)

% implemented - USPEX Version 8.5.0
global ORG_STRUC
global OFF_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numIons = ORG_STRUC.numIons;
[candidate, lat] = Random_Init_000(Ind_No, numIons);
OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';
%disp(['Structure ' num2str(Ind_No) ' generated with random symmetry: '...
%       num2str(nsym)]);
