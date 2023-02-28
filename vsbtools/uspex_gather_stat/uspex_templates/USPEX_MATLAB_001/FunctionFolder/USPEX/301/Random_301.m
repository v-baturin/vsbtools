function Random_301(Ind_No)

global ORG_STRUC
global OFF_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('Seeds/compositions')
    [nothing, result]=unix(['cat Seeds/compositions |wc -l']);
    tmp = str2num(result);
    if tmp(1)>0
        ORG_STRUC.firstGeneSplit=load('Seeds/compositions');
        ORG_STRUC.splitN = size(ORG_STRUC.firstGeneSplit,1);
    end
end

tmp = ceil(rand*ORG_STRUC.splitN);
if tmp == 0
    tmp = 1;
end
numBlocks = ORG_STRUC.firstGeneSplit(tmp,:);
numIons   = numBlocks*ORG_STRUC.numIons;
[candidate, lat] = Random_Init_301(Ind_No, numBlocks, 0);
OFF_STRUC.POPULATION(Ind_No).howCome = '  Random  ';
OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
OFF_STRUC.POPULATION(Ind_No).Parents = [];
OFF_STRUC.POPULATION(Ind_No).numIons = numIons;
OFF_STRUC.POPULATION(Ind_No).numBlocks = numBlocks;
