function [lat, new_Coord, deviation]= move_along_SoftMode_Mutation(coord, numIons, lattice, eigvector, mut_degree)
global ORG_STRUC

N = sum(numIons);
vec = zeros(1,3);
new_Coord = coord;
lat = lattice;

vect = zeros(1,3);
R_val = zeros(1,length(ORG_STRUC.atomType));
for i = 1 : length(ORG_STRUC.atomType)
 s = covalentRadius(ORG_STRUC.atomType(i));
 R_val(i) = str2num(s);
end
at_types = zeros(1,N);
for k = 1 : N
   tmp = k;
   while tmp > 0
     at_types(k) = at_types(k) + 1;
     tmp = tmp - numIons(at_types(k));
   end
end                                                  % END this part needed for clusters

if length(lattice) == 6
 lattice = latConverter(lattice);
end

coef = 0;
for i = 1 : N
  vec(1) = eigvector((i-1)*3+1);
  vec(2) = eigvector((i-1)*3+2);
  vec(3) = eigvector((i-1)*3+3);
  if norm(vec) > coef
   coef = norm(vec);
  end
end
normfac = mut_degree*ORG_STRUC.howManyMut/coef;   % this way the maximal displacement for an atom is equal to howManyMut*mut_degree
notDone = 1;
step = 0;

while notDone
 for i = 1 : N                      % shift the atoms
   vec(1) = eigvector((i-1)*3+1);
   vec(2) = eigvector((i-1)*3+2);
   vec(3) = eigvector((i-1)*3+3);
   new_Coord(i,:) = coord(i,:) + normfac*vec*inv(lattice);
   new_Coord(i,1) = new_Coord(i,1) - floor(new_Coord(i,1));
   new_Coord(i,2) = new_Coord(i,2) - floor(new_Coord(i,2));
   new_Coord(i,3) = new_Coord(i,3) - floor(new_Coord(i,3));
   deviation(i) = norm(normfac*vec);
 end
  notDone = 0;
 if step > 10
  notDone = 0;
 end
end
