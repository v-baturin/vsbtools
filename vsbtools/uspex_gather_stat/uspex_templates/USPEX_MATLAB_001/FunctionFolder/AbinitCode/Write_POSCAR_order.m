function Write_POSCAR_order(atomType, Count, symg, numIons, lattice, coor, order)



content = POSCARContent(atomType, Count, symg, numIons, lattice, coor, order);

writeContent2File('POSCAR_order', content, 'w');