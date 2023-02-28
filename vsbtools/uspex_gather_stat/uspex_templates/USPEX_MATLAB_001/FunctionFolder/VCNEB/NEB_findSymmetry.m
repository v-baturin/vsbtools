function NEB_findSymmetry()
%-------------------------------------------------%

global ORG_STRUC
global POP_STRUC
%-------------------------------------------------%


for i = 1:ORG_STRUC.numImages;
    lat = POP_STRUC.POPULATION(i).LATTICE;
    coord = POP_STRUC.POPULATION(i).COORDINATES;
    cd([ORG_STRUC.homePath '/CalcFoldTemp']);
    POP_STRUC.POPULATION(i).spacegroupNumber = anasym_stokes(lat, coord, ORG_STRUC.numIons, ORG_STRUC.atomType,ORG_STRUC.SGtolerance);
    cd([ORG_STRUC.homePath]);
end
