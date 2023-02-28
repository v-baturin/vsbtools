function Cell_transformation_201(whichInd)
global ORG_STRUC
global POP_STRUC

Step     = POP_STRUC.POPULATION(whichInd).Step;
N_Step   = length([ORG_STRUC.abinitioCode]);

lattice    = POP_STRUC.POPULATION(whichInd).LATTICE;
coordinate = POP_STRUC.POPULATION(whichInd).COORDINATES;
numIons    = POP_STRUC.POPULATION(whichInd).numIons;

[lat,candidate,numIons]=reduce_Surface(lattice,coordinate,numIons,whichInd,1);
%disp('surfaceafter GULP')
%disp('1.0000')
%disp(num2str(lattice))
%disp('Al O')
%disp(num2str(numIons))
%disp('direct')
%disp(num2str(candidate))
if size(candidate,1)==0
    INIT_numIons = POP_STRUC.POPULATION(whichInd).INIT_numIons;
    Bulk_numIons = POP_STRUC.POPULATION(whichInd).Bulk_numIons;
    if sum(INIT_numIons - Bulk_numIons) > 0 % if the initial structure is not bare surface
       disp('No add atoms now, start to dock')
       GenerateSurface(whichInd, 1);
    end
else
    POP_STRUC.POPULATION(whichInd).Surface_LATTICE=lat;
    POP_STRUC.POPULATION(whichInd).Surface_COORDINATES=candidate;
    POP_STRUC.POPULATION(whichInd).Surface_numIons=numIons;
end

if Step < N_Step
   surcoor    = POP_STRUC.POPULATION(whichInd).Surface_COORDINATES;
   surlat     = POP_STRUC.POPULATION(whichInd).Surface_LATTICE;
   surnumIons = POP_STRUC.POPULATION(whichInd).Surface_numIons;
   bulklat    = POP_STRUC.POPULATION(whichInd).Bulk_LATTICE;
   bulkcoor   = POP_STRUC.POPULATION(whichInd).Bulk_COORDINATES;
   ntyp=POP_STRUC.POPULATION(whichInd).Bulk_numIons;
   [lat,candidate,numIons,chanAList] = makeSurface(surlat,surcoor,surnumIons,...
                         bulklat,bulkcoor,ntyp,ORG_STRUC.vacuumSize(Step));
   POP_STRUC.POPULATION(whichInd).numIons     = numIons;
   POP_STRUC.POPULATION(whichInd).LATTICE     = lat;
   POP_STRUC.POPULATION(whichInd).COORDINATES = candidate;
   POP_STRUC.POPULATION(whichInd).chanAList   = chanAList;
end
