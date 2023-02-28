function WriteGenerationOutput_311(fitness)

% USPEX Version 9.3.2
% Change: new output, hardness, order, etc

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC
global POOL_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resFolder = POP_STRUC.resFolder;
convex_hull = POP_STRUC.convex_hull;
atomType    = ORG_STRUC.atomType;

IND = [];
if ~isempty(convex_hull)
    fpath1 = [ resFolder '/convex_hull'];
    fp1 = fopen(fpath1, 'a+');
    fprintf(fp1, '---- generation%3d ----\n', POP_STRUC.generation);
    item = 1;
    for i = 1:size(convex_hull,1)
        IND(i) = POP_STRUC.convex_hull(i,end);
        
        for j=1:length(atomType)
            fprintf(fp1,'%4d', POP_STRUC.POPULATION(IND(i)).numIons(j));
        end
        fprintf(fp1,'%12.4f\n', convex_hull(i,end-1));
    end
    fclose(fp1);
else
    IND=POP_STRUC.ranking(1);
end

for i=1:length(IND)
    lattice = POP_STRUC.POPULATION(IND(i)).LATTICE;
    coor    = POP_STRUC.POPULATION(IND(i)).COORDINATES;
    numIons = POP_STRUC.POPULATION(IND(i)).numIons;
    lattice = POP_STRUC.POPULATION(IND(i)).LATTICE;
    symg    = POP_STRUC.POPULATION(IND(i)).symg;
    order   = POP_STRUC.POPULATION(IND(i)).order;
    count   = POP_STRUC.POPULATION(IND(i)).Number;
    Write_POSCAR(atomType, count, symg, numIons, lattice, coor);
    Write_POSCAR_order(atomType, count, symg, numIons, lattice, coor, order);
    
    unixCmd([' cat POSCAR       >>' resFolder '/BESTgatheredPOSCARS']);
    unixCmd([' cat POSCAR_order >>' resFolder '/BESTgatheredPOSCARS_order']);
end

update_USPEX_GENERATION(IND, fitness, 2);

% FIGURES
cd (resFolder)
if ~isempty(POP_STRUC.convex_hull)
    extendedConvexHull_311(POP_STRUC.convex_hull, 0.5);
    if ORG_STRUC.PhaseDiagram
        phase_diagram_varcomp;
    end
else
    extractStructures(ORG_STRUC.populationSize, ORG_STRUC.numIons);
end
makeFigures(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
%update VOLUME
V = [];
comp = [];
for i=1:length(USPEX_STRUC.POPULATION)
    V = [V; det(USPEX_STRUC.POPULATION(i).LATTICE)];
    comp = [comp; USPEX_STRUC.POPULATION(i).numBlocks];
end
ORG_STRUC.latVolume = (comp\V)';
unixCmd(['echo "Approximate volume(s): "' num2str(ORG_STRUC.latVolume) ' >> OUTPUT.txt']);

cd ..

% update compositions
comp = load('Seeds/compositions');
fp = fopen('Seeds/compositions', 'w');
for i = 1:size(comp,1)
    towrite = 1;
    block = comp(i,:)*ORG_STRUC.numIons;
    ratio = block/sum(block);
    for j=1:length(POOL_STRUC.Composition_ranking)
        if sum(abs(ratio-POOL_STRUC.Composition_ratio(j,:))) < 0.01
            if POOL_STRUC.Composition_ranking(j)<0
                towrite = 0;
            end
            break;
        end
    end
    if towrite==1
        if CompositionCheck(block) %from AntiSeeds
            for k=1:size(comp,2)
                fprintf(fp, '%4d', comp(i,k));
            end
            fprintf(fp, '\n');
        end
    end
end
fclose(fp);

WriteIndividual(resFolder);
WriteBest(resFolder);
WriteGeneration(resFolder);
WriteProperties(ORG_STRUC.optType, resFolder);
WriteCompStatistic(resFolder)

safesave([resFolder '/POOL.mat'], POOL_STRUC);
safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
