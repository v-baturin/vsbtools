function ID = update_USPEX_GENERATION(IND, varcomp, Step)

%To update all the necessary items at the end of each gen to USPEX_STRUC
%Output the BESTINDIVIDUALS info
%Lastly updated by Qiang Zhu (2014/02/18)

global USPEX_STRUC
global POP_STRUC
resFolder = POP_STRUC.resFolder;

generation = POP_STRUC.generation + 1;

ID = POP_STRUC.POPULATION(IND).Number;
N_atom = sum(POP_STRUC.POPULATION(IND).numIons);
USPEX_STRUC.GENERATION(generation).BestID = ID;
USPEX_STRUC.GENERATION(generation).Best_enth = POP_STRUC.POPULATION(IND).Enthalpies(Step)/N_atom;
USPEX_STRUC.GENERATION(generation).Best_enth_relaxed = POP_STRUC.POPULATION(IND).Enthalpies(end)/N_atom;
%%%%%%lattice_basic.dat
lat_basic = POP_STRUC.lat; 
lat_basic(1,:) = lat_basic(1,:)/POP_STRUC.POPULATION(IND).superCell(1);
lat_basic(2,:) = lat_basic(2,:)/POP_STRUC.POPULATION(IND).superCell(2);
lat_basic(3,:) = lat_basic(3,:)/POP_STRUC.POPULATION(IND).superCell(3);
lat_basic1 = MattoVec(lat_basic);
USPEX_STRUC.GENERATION(generation).lat_basic = lat_basic1;

unixCmd(['echo ' num2str(lat_basic1,'%12.4f ') '  >> ' resFolder '/lattice_basic.dat']);

if varcomp
   if sum(POP_STRUC.POPULATION(IND).superCell) == 3 % [1 1 1]
      POP_STRUC.bestBasicStructure = POP_STRUC.POPULATION(IND).COORDINATES;
      USPEX_STRUC.GENERATION(generation).bestBasicStructure = POP_STRUC.bestBasicStructure;
   end
end
