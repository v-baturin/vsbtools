function RandTop_301(Ind_No)

global ORG_STRUC
global OFF_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using topology %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if exist('Seeds/compositions')
   %[a, b]=unix(['wc Seeds/compositions']); %Sometimes the composition file is broken
   %tmp = str2num(b(1:2));
   [nothing, result]=unix(['cat Seeds/compositions |wc -l']);
   tmp = str2num(result);
   if tmp(1)>0
   ORG_STRUC.firstGeneSplit=load('Seeds/compositions');
   ORG_STRUC.splitN = length(ORG_STRUC.firstGeneSplit);
   end
end   

tmp = ceil(rand*ORG_STRUC.splitN);
if tmp == 0
    tmp = 1;
end
numBlocks = ORG_STRUC.firstGeneSplit(tmp,:);

[candidate, lat] = Random_Init_301(Ind_No, numBlocks,1);%1 - topological random

OFF_STRUC.POPULATION(Ind_No).COORDINATES = candidate;
OFF_STRUC.POPULATION(Ind_No).LATTICE = lat;
OFF_STRUC.POPULATION(Ind_No).Parents = [];
OFF_STRUC.POPULATION(Ind_No).numIons = numBlocks*ORG_STRUC.numIons;
OFF_STRUC.POPULATION(Ind_No).numBlocks = numBlocks;
OFF_STRUC.POPULATION(Ind_No).howCome = '  RandTop  ';

