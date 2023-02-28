function WriteIndividualOutput_M400(Ind_No)
% $Rev$
% $Author$
% $Date$

global POP_STRUC
global ORG_STRUC

resFolder = POP_STRUC.resFolder;
Count=POP_STRUC.POPULATION(Ind_No).Number;

if (Ind_No == 0)
    bodyCount = 0;
else
    bodyCount = POP_STRUC.POPULATION(Ind_No).Number;
end

if strcmp(POP_STRUC.POPULATION(Ind_No).howCome, 'keptBest') > 0
    cd(resFolder);

    parent_id = str2num(POP_STRUC.POPULATION(Ind_No).Parents.parent);

    unixCmd(['cat /dev/null > PDB']);
    getPDB(parent_id);
    unixCmd(['cat PDB | sed "s/HEADER    EA' num2str(parent_id) '/HEADER    EA' num2str(bodyCount) ' <-' num2str(parent_id) '/g" >> gatheredPDB']);
    unixCmd(['rm -f PDB']);

    unixCmd(['cat /dev/null > POSCAR']);
    getPOSCAR(parent_id, 'gatheredPOSCARS');
    unixCmd(['cat POSCAR | sed "s/EA' num2str(parent_id) '/EA' num2str(bodyCount) '/g" >> gatheredPOSCARS']);
    unixCmd(['rm -f POSCAR']);

    cd ..
else
    unixCmd(['cat PDB | sed "s/HEADER    /HEADER    EA' num2str(bodyCount) ' /g" >> ' POP_STRUC.resFolder '/gatheredPDB']);

    unixCmd(['echo EA' num2str(bodyCount) ' >> ' POP_STRUC.resFolder '/gatheredMAKE']);
    unixCmd(['cat MAKE                      >> ' POP_STRUC.resFolder '/gatheredMAKE']);

    unixCmd(['cat POSCAR     | sed "s/EA0000/EA' num2str(bodyCount) '/g" >> ' POP_STRUC.resFolder '/gatheredPOSCARS']);
end

unixCmd(['rm -f PDB MAKE POSCAR']);

atomType = ORG_STRUC.atomType;
update_USPEX_INDIVIDUAL(POP_STRUC.POPULATION(Ind_No), resFolder, ...
                        POP_STRUC.generation, atomType);

WriteOUTPUT(Count, resFolder);

%INIT_numIons=POP_STRUC.POPULATION(Ind_No).INIT_numIons;
%INIT_LAT    =POP_STRUC.POPULATION(Ind_No).INIT_LAT;
%INIT_COORD  =POP_STRUC.POPULATION(Ind_No).INIT_COORD;
%Write_InitStructure(Count, INIT_numIons, INIT_LAT, INIT_COORD, resFolder);
