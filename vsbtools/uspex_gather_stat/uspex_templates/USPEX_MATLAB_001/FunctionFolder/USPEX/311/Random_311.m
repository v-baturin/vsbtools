function Random_311(Ind_No)

% implemented - USPEX Version 9.3.7
global ORG_STRUC
global OFF_STRUC

nSplit = length( ORG_STRUC.firstGeneSplit );
numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) , :);

k=1;
while ~CompositionCheck( numBlocks )
    numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) , :);
    k = k+1;
    if k > 1000
        disp('Error in Anti-compositions file, please check !!');
        quit
    end
end

numMols = numBlocks*ORG_STRUC.numIons;
POP     = Random_Init_311(Ind_No, numMols);

OFF_STRUC.POPULATION(Ind_No).MOLECULES  =  POP.MOLECULES;
OFF_STRUC.POPULATION(Ind_No).numMols    =  POP.numMols;
OFF_STRUC.POPULATION(Ind_No).MtypeLIST  =  POP.MtypeLIST;
OFF_STRUC.POPULATION(Ind_No).typesAList =  POP.typesAList;
OFF_STRUC.POPULATION(Ind_No).numIons    =  POP.numIons;
OFF_STRUC.POPULATION(Ind_No).LATTICE    =  POP.LATTICE;
OFF_STRUC.POPULATION(Ind_No).howCome    =  ' Random ';
OFF_STRUC.POPULATION(Ind_No).numBlocks  = numBlocks;
