function WriteGenerationOutput_300(fitness)

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atomType  = ORG_STRUC.atomType;
resFolder = POP_STRUC.resFolder;


if ORG_STRUC.paretoRanking ~= 0
    IND=POP_STRUC.ranking(1:POP_STRUC.paretoFront(1));
else
    IND=POP_STRUC.ranking(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Make Sure bodyCount is correct
for i=1:length(IND)
    lattice = POP_STRUC.POPULATION(IND(i)).LATTICE;
    coor    = POP_STRUC.POPULATION(IND(i)).COORDINATES;
    numIons = POP_STRUC.POPULATION(IND(i)).numIons;
    symg    = POP_STRUC.POPULATION(IND(i)).symg;
    order   = POP_STRUC.POPULATION(IND(i)).order;
    count   = POP_STRUC.POPULATION(IND(i)).Number;
    Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
    Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);

    unixCmd([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
    unixCmd([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);
end
update_USPEX_GENERATION(IND, fitness, 1);
% FIGURES
cd (resFolder);
makeFigures(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
extractStructures(ORG_STRUC.populationSize, ORG_STRUC.weight);
if ORG_STRUC.PhaseDiagram
    phase_diagram;
end

V = [];
for i=1:length(USPEX_STRUC.POPULATION)
    V = [V det(USPEX_STRUC.POPULATION(i).LATTICE)*sum(ORG_STRUC.numIons)/sum(USPEX_STRUC.POPULATION(i).numIons) ];
end
ORG_STRUC.latVolume = mean(V);

cd ..
fpath = [ resFolder '/OUTPUT.txt'];
fp = fopen(fpath, 'a+');
fprintf(fp, [alignLine( sprintf('Approximate volume(s): %.4f A^3', mean(V)) ) '\n'] );
fclose(fp);

WriteIndividual(resFolder);
WriteBest(resFolder);
WriteGeneration(resFolder);
WriteProperties(ORG_STRUC.optType, resFolder);
WriteTEproperties(resFolder);
WriteCompStatistic(resFolder);
if ORG_STRUC.paretoRanking ~= 0
    WriteParetoRanking(resFolder);
end

safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
