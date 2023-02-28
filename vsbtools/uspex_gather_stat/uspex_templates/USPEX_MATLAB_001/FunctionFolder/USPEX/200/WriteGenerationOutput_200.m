function WriteGenerationOutput_200(fitness)

% USPEX Version 9.3.2
% Change: new output, hardness, order, etc

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atomType= ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;

IND=POP_STRUC.ranking(1);
lattice = POP_STRUC.POPULATION(IND).LATTICE;
coor    = POP_STRUC.POPULATION(IND).COORDINATES;
numIons = POP_STRUC.POPULATION(IND).numIons;
lattice = POP_STRUC.POPULATION(IND).LATTICE;
symg    = POP_STRUC.POPULATION(IND).symg;
%order   = POP_STRUC.POPULATION(IND).order;
count   = POP_STRUC.POPULATION(IND).Number;

count1 = 0;
order = zeros(sum(numIons),1);
Surface_order   = POP_STRUC.POPULATION(IND).Surface_order;
for i=1:size(coor,1)
    if POP_STRUC.POPULATION(IND).chanAList(i)==1
       count1=count1+1;
       order(i) = Surface_order(count1);
    end
end

Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
unixCmd([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
unixCmd([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);
update_USPEX_GENERATION(IND, fitness, 1);

% FIGURES
cd(resFolder);
makeFigures(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
extractStructures(20, ORG_STRUC.weight);
cd ..
WriteIndividual(resFolder);
WriteBest(resFolder);
WriteGeneration(resFolder);

if exist('writeMagmoment')
    writeMagmoment(resFolder);
end

safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
