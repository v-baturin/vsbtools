function Cell_transformation_M200(whichInd)
global ORG_STRUC
global POP_STRUC

Step     = POP_STRUC.POPULATION(whichInd).Step;
N_Step   = length([ORG_STRUC.abinitioCode]);

lattice    = POP_STRUC.POPULATION(whichInd).LATTICE;
coordinate = POP_STRUC.POPULATION(whichInd).COORDINATES;
[lat,candidate]=reduce2D(lattice,coordinate,ORG_STRUC.thicknessS+3);

POP_STRUC.POPULATION(whichInd).LATTICE_2D     = lat;
POP_STRUC.POPULATION(whichInd).COORDINATES_2D = candidate;
POP_STRUC.POPULATION(whichInd).Volume_2D      = det(lat);

if Step < N_Step
   [lattice, coordinate] = make2D(lat, candidate, ORG_STRUC.vacuumSize(Step)-3);
   POP_STRUC.POPULATION(whichInd).LATTICE = lattice;
   POP_STRUC.POPULATION(whichInd).COORDINATES = coordinate;
else
   if max(coordinate(:,3))-min(coordinate(:,3)) > (1-2/lattice(3,3)) %too thick
      POP_STRUC.POPULATION(whichInd).Enthalpies(Step-1)=100000;
   end
end
