function Cell_transformation_200(whichInd)
global ORG_STRUC
global POP_STRUC

Step     = POP_STRUC.POPULATION(whichInd).Step;
N_Step   = length([ORG_STRUC.abinitioCode]);

lattice    = POP_STRUC.POPULATION(whichInd).LATTICE;
coordinate = POP_STRUC.POPULATION(whichInd).COORDINATES;
numIons    = POP_STRUC.POPULATION(whichInd).numIons;
cell       = POP_STRUC.POPULATION(whichInd).cell;
[lat,candidate,numIons]=reduce_Surface(lattice,coordinate,numIons,whichInd,0);
if isequal(numIons, ORG_STRUC.numIons*det(reshape(cell,[2,2])) )
   POP_STRUC.POPULATION(whichInd).Surface_LATTICE=lat;
   POP_STRUC.POPULATION(whichInd).Surface_COORDINATES=candidate;
   POP_STRUC.POPULATION(whichInd).Surface_numIons=numIons;
elseif sum(numIons)==0
   GenerateSurface(whichInd, 0);
else
   MakeupSurface(whichInd,lat,candidate,numIons);
end

if Step < N_Step %At final step, we don't make new surface
   surlat     = POP_STRUC.POPULATION(whichInd).Surface_LATTICE;
   surcoor    = POP_STRUC.POPULATION(whichInd).Surface_COORDINATES;
   surnumIons = POP_STRUC.POPULATION(whichInd).Surface_numIons;
   bulklat    = POP_STRUC.POPULATION(whichInd).Bulk_LATTICE;
   bulkcoor   = POP_STRUC.POPULATION(whichInd).Bulk_COORDINATES;
   ntyp       = POP_STRUC.POPULATION(whichInd).Bulk_numIons;
   [lat,candidate,numIons,chanAList] =  makeSurface...
   (surlat,surcoor,surnumIons,bulklat,bulkcoor,ntyp,ORG_STRUC.vacuumSize(Step));
   POP_STRUC.POPULATION(whichInd).numIons = numIons;
   POP_STRUC.POPULATION(whichInd).LATTICE = lat;
   POP_STRUC.POPULATION(whichInd).COORDINATES = candidate;
   POP_STRUC.POPULATION(whichInd).chanAList=chanAList;
end
