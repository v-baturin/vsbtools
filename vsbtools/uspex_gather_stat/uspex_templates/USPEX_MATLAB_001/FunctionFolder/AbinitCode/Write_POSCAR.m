function Write_POSCAR(atomType, Count, symg, numIons, lattice, coor)


content = POSCARContent(atomType, Count, symg, numIons, lattice, coor);

writeContent2File('POSCAR', content, 'w');
%        if POP_STRUC.POPULATION(whichInd).Step > ORG_STRUC.conv_till && ORG_STRUC.conv_till < length([ORG_STRUC.abinitioCode])
%        copyfile('symmetrizedPOSCAR', 'POSCAR')
