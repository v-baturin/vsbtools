function STM4_WriteOutFile(opt)

global POP_STRUC
global ORG_STRUC

%--------------------------------------------------------------------%
atomType  = ORG_STRUC.atomType;
numImages = ORG_STRUC.numImages;
numIons  = ORG_STRUC.numIons;


if     0==opt
    file = [ ORG_STRUC.resFolder '/PATH/path_' num2str(POP_STRUC.step) '.POSCAR'];
    opt='w';
elseif 1==opt
    file = [ ORG_STRUC.resFolder '/transitionPath_POSCARs' ];
    opt='w';
elseif 2==opt
    file = [ ORG_STRUC.resFolder, '/AuxiliaryFiles/gatheredPOSCARS'];;
    opt='a';
end

%--------------------------------------------------------------------------------
content={};
for i = 1:numImages
    if 2==opt
        header = sprintf('Step %3d -- Image %-3d', POP_STRUC.step, i);
    else
        header = sprintf('Image %-3d', i);
    end

    newCont = {};
    lattice = POP_STRUC.POPULATION(i).LATTICE;
    coord   = POP_STRUC.POPULATION(i).COORDINATES;
 
    newCont = POSCARContent(atomType, header, '', numIons, lattice, coord);
    if i==1
       content = newCont;
    else
       content(end+1:end+length(newCont)) = newCont; 
    end
end

writeContent2File(file, content, opt);

