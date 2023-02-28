function WriteGenerationOutput_PSO(fitness)

% USPEX Version 9.3.2
% Change: new output, hardness, order, etc

global POP_STRUC
global ORG_STRUC
global PSO_STRUC
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

Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
unixCmd([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
unixCmd([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);

update_PSO_GENERATION(IND, fitness, 1);

% FIGURES
cd (resFolder);
makeFigures(ORG_STRUC.pickUpGen, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
extractStructures_PSO(ORG_STRUC.populationSize, ORG_STRUC.weight);
cd ..
V = [];
for i=1:length(PSO_STRUC.POPULATION)
    V = [V det(PSO_STRUC.POPULATION(i).LATTICE)];
end
ORG_STRUC.latVolume = mean(V);
WriteIndividual(resFolder);
WriteBest(resFolder);
WriteGeneration(resFolder);
WriteProperties(ORG_STRUC.optType, resFolder);
safesave([resFolder '/PSO.mat'], PSO_STRUC);

