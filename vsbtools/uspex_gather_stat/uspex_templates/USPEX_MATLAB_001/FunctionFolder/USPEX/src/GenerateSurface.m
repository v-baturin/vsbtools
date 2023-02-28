function GenerateSurface(whichInd, varcomp)

global ORG_STRUC
global POP_STRUC 

CellList = load('Seeds/cell');

Vacuum = ORG_STRUC.vacuumSize(1);
ID = ceil(rand()*size(CellList,1));
cell=CellList(ID,:);
if varcomp
   ANS  = Random_Init_201(whichInd, Vacuum, cell);
else
   ANS  = Random_Init_200(whichInd, Vacuum, cell);
end

POP_STRUC.POPULATION(whichInd).COORDINATES         = ANS(1).candidate;
POP_STRUC.POPULATION(whichInd).numIons             = ANS(1).numIons;
POP_STRUC.POPULATION(whichInd).LATTICE             = ANS(1).lat;
POP_STRUC.POPULATION(whichInd).chanAList           = ANS(1).chanAList;
POP_STRUC.POPULATION(whichInd).Surface_LATTICE     = ANS(1).sur_lat;
POP_STRUC.POPULATION(whichInd).Surface_COORDINATES = ANS(1).sur_candidate;
POP_STRUC.POPULATION(whichInd).Surface_numIons     = ANS(1).sur_numIons;
POP_STRUC.POPULATION(whichInd).Bulk_LATTICE        = ANS(1).bulk_lat;
POP_STRUC.POPULATION(whichInd).Bulk_COORDINATES    = ANS(1).bulk_pos;
POP_STRUC.POPULATION(whichInd).Bulk_numIons        = ANS(1).bulk_numIons;
POP_STRUC.POPULATION(whichInd).cell                = cell;
POP_STRUC.POPULATION(whichInd).howCome             = '  Random  ';
POP_STRUC.POPULATION(whichInd).Step                = 1;

