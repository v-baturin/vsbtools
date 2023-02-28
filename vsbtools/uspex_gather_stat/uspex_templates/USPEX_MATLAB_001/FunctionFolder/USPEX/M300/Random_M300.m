function Random_M300(Ind_No)

% implemented - USPEX Version 8.5.0

global ORG_STRUC
global OFF_STRUC

nSplit = length( ORG_STRUC.firstGeneSplit );
numBlocks = ORG_STRUC.firstGeneSplit( ceil(rand(1)*nSplit) );
numIons   = ORG_STRUC.numIons*numBlocks;
POP = Random_Init_M300(Ind_No, numIons);
OFF_STRUC.POPULATION(Ind_No).Bulk_LATTICE     = POP.Bulk_LATTICE;
OFF_STRUC.POPULATION(Ind_No).Bulk_COORDINATES = POP.Bulk_COORDINATES;
OFF_STRUC.POPULATION(Ind_No).Bulk_typesAList  = POP.Bulk_typesAList;
OFF_STRUC.POPULATION(Ind_No).Bulk_numIons     = POP.Bulk_numIons;
OFF_STRUC.POPULATION(Ind_No).COORDINATES      = POP.COORDINATES;
OFF_STRUC.POPULATION(Ind_No).numIons          = POP.numIons;         
OFF_STRUC.POPULATION(Ind_No).LATTICE          = POP.LATTICE; 
OFF_STRUC.POPULATION(Ind_No).typesAList       = POP.typesAList;
OFF_STRUC.POPULATION(Ind_No).GB_LATTICE       = POP.GB_LATTICE;      
OFF_STRUC.POPULATION(Ind_No).GB_COORDINATES   = POP.GB_COORDINATES;
OFF_STRUC.POPULATION(Ind_No).GB_numIons       = POP.GB_numIons;
OFF_STRUC.POPULATION(Ind_No).GB_typesAList    = POP.GB_typesAList;  
OFF_STRUC.POPULATION(Ind_No).chanAList        = POP.chanAList;
OFF_STRUC.POPULATION(Ind_No).howCome          = ' Random ';
