function WriteGenerationOutput_M400(fitness)
% $Rev$
% $Author$
% $Date$

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% OUTPUT THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resFolder = POP_STRUC.resFolder;

IND=POP_STRUC.ranking(1);

update_USPEX_GENERATION(IND, fitness, 1);

cd(resFolder);

getPDB(POP_STRUC.POPULATION(IND).Number);
unixCmd(['cat PDB >> BESTgatheredPDB']);

pdb_folder = 'BEST_PDB';
if ~isequal(exist(pdb_folder, 'dir'), 7) % 7 = directory.
    mkdir(pdb_folder);
end
unixCmd(['cp -f PDB ' pdb_folder '/EA' sprintf('%04d', POP_STRUC.POPULATION(IND).Number) '.pdb']);
unixCmd(['rm -f PDB']);

getPOSCAR(POP_STRUC.POPULATION(IND).Number, 'gatheredPOSCARS');
unixCmd(['cat POSCAR >> BESTgatheredPOSCARS']);
unixCmd(['rm -f POSCAR']);

% FIGURES
makeFigures_M400(ORG_STRUC.pickUpNCount, length(ORG_STRUC.abinitioCode), ORG_STRUC.constLattice);
cd ..

WriteIndividual(resFolder, 'kcal/mol');
WriteBest(resFolder, 'kcal/mol');
WriteGeneration(resFolder);

safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
