function MINHOP_WriteGenerationOutput(fitness, flag)

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%flag = 1:  if there is a FullRelax tag, and we need to keep the BESTgatheredPOSCARS
%flag = 0:  OUTPUT everything
atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;
IND=POP_STRUC.ranking(1);

lattice = POP_STRUC.POPULATION(IND).LATTICE;
coor    = POP_STRUC.POPULATION(IND).COORDINATES;
numIons = POP_STRUC.POPULATION(IND).numIons;
lattice = POP_STRUC.POPULATION(IND).LATTICE;
symg    = POP_STRUC.POPULATION(IND).symg;
order   = POP_STRUC.POPULATION(IND).order;
count   = POP_STRUC.POPULATION(IND).Number;

%%%%%%POSCARs
if flag == 1
   Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
   Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
   unixCmd([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
   unixCmd([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);

else
   %%%%%%BESTIndividuals
   ID = update_USPEX_GENERATION(IND, ORG_STRUC.varcomp, ORG_STRUC.conv_till);
   MINHOP_WriteOUTPUT(ID, resFolder, ORG_STRUC.varcomp, ORG_STRUC.conv_till, 3);
       Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
       unixCmd([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS_relaxed']);
       MINHOP_WriteOUTPUT(ID, resFolder, ORG_STRUC.varcomp, length([ORG_STRUC.abinitioCode]), 4);
   %%%%%%
   safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
   %% FIGURES %%
   cd (resFolder);
%   MINHOP_makeFigures(ORG_STRUC.conv_till);
   if (ORG_STRUC.FullRelax > 0)
   MINHOP_extractStructures(round(ORG_STRUC.numGenerations/10), ORG_STRUC.weight);
   end
   cd ..
end
