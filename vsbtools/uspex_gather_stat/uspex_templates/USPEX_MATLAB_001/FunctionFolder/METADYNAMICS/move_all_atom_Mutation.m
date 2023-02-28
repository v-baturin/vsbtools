function [new_Coord] = move_all_atom_Mutation(coord, numIons, lattice, max_sigma)

N = sum(numIons);
new_Coord = coord;

if length(lattice) == 6
 lattice = latConverter(lattice);
end
temp_potLat = latConverter(lattice);

 for i = 1 : N
  deviat_dist = randn(3,1)*max_sigma;

  new_Coord(i,1) = coord(i,1) + deviat_dist(1)/temp_potLat(1);
  new_Coord(i,2) = coord(i,2) + deviat_dist(2)/temp_potLat(2);
  new_Coord(i,3) = coord(i,3) + deviat_dist(3)/temp_potLat(3);

  new_Coord(i,1) = new_Coord(i,1) - floor(new_Coord(i,1));
  new_Coord(i,2) = new_Coord(i,2) - floor(new_Coord(i,2));
  new_Coord(i,3) = new_Coord(i,3) - floor(new_Coord(i,3));
 end
