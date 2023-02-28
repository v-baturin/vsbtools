function WriteIndividualOutput_PSO(Ind_No)
% USPEX 9.4.0
% We choose to update the files requiring I/O files first
% In case sometimes I/O errors in file system
% Lastly updated by Qiang Zhu (2013/11/22)

global POP_STRUC
global ORG_STRUC

atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;

symg        =POP_STRUC.POPULATION(Ind_No).symg;
count       =POP_STRUC.POPULATION(Ind_No).Number;
INIT_numIons=POP_STRUC.POPULATION(Ind_No).INIT_numIons;
INIT_LAT    =POP_STRUC.POPULATION(Ind_No).INIT_LAT;
INIT_COORD  =POP_STRUC.POPULATION(Ind_No).INIT_COORD;
Write_POSCAR(atomType, count, symg, INIT_numIons, INIT_LAT, INIT_COORD);
unixCmd([' cat POSCAR      >> ' POP_STRUC.resFolder '/gatheredPOSCARS_unrelaxed']);


lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
coor    = POP_STRUC.POPULATION(Ind_No).COORDINATES;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
order   = POP_STRUC.POPULATION(Ind_No).order;

% write symmetrized structure in .cif format
if ORG_STRUC.doSpaceGroup == 1
    current_path = pwd;
    cd([ORG_STRUC.homePath '/CalcFoldTemp']);
    POP_STRUC.POPULATION(Ind_No).symg = anasym_stokes(lattice, coor, numIons, atomType, ORG_STRUC.SGtolerance);
    symg        =POP_STRUC.POPULATION(Ind_No).symg;
    cd(current_path)
    unixCmd(['echo data_findsym-STRUC-' num2str(count) '            >> ' resFolder '/symmetrized_structures.cif']);
    unixCmd(['cat CalcFoldTemp/symmetrized.cif                    >> ' resFolder '/symmetrized_structures.cif']);
end


Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
unixCmd([' cat POSCAR      >> ' POP_STRUC.resFolder '/gatheredPOSCARS']);
unixCmd([' cat POSCAR_order >>' POP_STRUC.resFolder '/gatheredPOSCARS_order']);


update_PSO_INDIVIDUAL(POP_STRUC.POPULATION(Ind_No), resFolder, ...
    POP_STRUC.generation, ORG_STRUC.atomType);

WriteOUTPUT_PSO(count, resFolder);

