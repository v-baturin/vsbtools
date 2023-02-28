function [lat1, coord1, order1, numIons1]=cresupercell(lattice, coord, numIons, lat_base, cell, atomType)

lat1 = [];
coord1 = [];
order1 = [];
numIons1 = [];
if isempty(coord)
    return
else 
   %1, Pick 1*1 cell
   [lat,coor,numIons] = Pick_primitive(lattice, coord, numIons, lat_base, atomType);
   
   %2, 
   if ~isempty(coor)
      %disp('cresupercell')
      %disp(['cresupercell: numIons: ' num2str(numIons)])
      [coor, numIons] = RemoveDup(coor, numIons, 0.05);
      %disp(['cresupercell: numIons: ' num2str(numIons)])
      [coord1, lat1, numIons1, Trans] = Make_MultiCell(coor, lat, numIons, cell);
      order1 = ones(sum(numIons1),1);
      %[Ni, V, dist_matrix, typ_i, typ_j, ho, ht] = makeMatrices_2D...
      %                               (lat1, coord1, sum(numIons1), 1);
      %[order1, FINGERPRINT, atom_fing] = fingerprint_calc_2D...
      %    (Ni, V, dist_matrix, typ_i, typ_j, sum(numIons1), ho, ht);
   end
end
