function Start_POP()

global ORG_STRUC
global POP_STRUC

N_Step = length([ORG_STRUC.abinitioCode]);
POP_STRUC.DoneOrder = zeros(1, length(POP_STRUC.POPULATION));
% when = 1, we do 'final' property calculation Haiyang Niu 2-17-2016
POP_STRUC.finalOptimization = 0; 
fractions = {'fracRand', 'fracRandTop', 'fracGene', 'fracG1ene', ...
             'fracPerm', 'fracAddAtom', 'fracRemAtom',...
             'fracTrans','fracLatMut',...
             'fracAtomsMut','fracRotMut','fracSpin',...
             'fracSecSwitch', 'fracShiftBorder'};  

for i = 1:length(fractions)
    eval( ['POP_STRUC.' fractions{i} '= ORG_STRUC.' fractions{i} ';'] );
end

for i = 1:length(POP_STRUC.POPULATION)
    if isempty(POP_STRUC.POPULATION(i).Step)
       POP_STRUC.POPULATION(i).Step = 1;
       POP_STRUC.POPULATION(i).Enthalpies = 100000*ones(1, N_Step);
       POP_STRUC.POPULATION(i).K_POINTS = ones(N_Step, 3);
       if ORG_STRUC.opt_sign == -1
          POP_STRUC.POPULATION(i).dielectric_tensor =100000*ones(1,6);
          POP_STRUC.POPULATION(i).gap = 100000;
          POP_STRUC.POPULATION(i).mag_moment = 0;
          POP_STRUC.POPULATION(i).TE_property = 100000;
	  POP_STRUC.POPULATION(i).Fphon = 0;
       else
          POP_STRUC.POPULATION(i).dielectric_tensor =zeros(1,6);
          POP_STRUC.POPULATION(i).gap = 0;
          POP_STRUC.POPULATION(i).hardness = 0;
          POP_STRUC.POPULATION(i).mag_moment = 100000;
          POP_STRUC.POPULATION(i).TE_property = 0;
	  POP_STRUC.POPULATION(i).Fphon = 100000;
       end
    end
    POP_STRUC.POPULATION(i).Error = 0;
    POP_STRUC.POPULATION(i).Folder = 0;
    POP_STRUC.POPULATION(i).ToDo = 1;
    POP_STRUC.POPULATION(i).Done = 0;
    POP_STRUC.POPULATION(i).JobID = 0;
    POP_STRUC.POPULATION(i).Number = 0;
    if ORG_STRUC.molecule == 1 
       POP_STRUC.POPULATION(i).COORDINATES = Sort_Atom(ORG_STRUC.atomType,...
                                        POP_STRUC.POPULATION(i).MOLECULES,...
                                        POP_STRUC.POPULATION(i).LATTICE,  ...
                                        POP_STRUC.POPULATION(i).typesAList);
    end
    POP_STRUC.POPULATION(i).INIT_COORD = POP_STRUC.POPULATION(i).COORDINATES;
    POP_STRUC.POPULATION(i).INIT_LAT   = POP_STRUC.POPULATION(i).LATTICE;
    POP_STRUC.POPULATION(i).INIT_numIons=POP_STRUC.POPULATION(i).numIons;
end

ReRank();
for indic = 1:ORG_STRUC.numParallelCalcs
    POP_STRUC.CalcFold(indic) = indic;
end
POP_STRUC.CalcFold_max = indic;

function COORDINATES = Sort_Atom(atomType, MOLECULES, LATTICE, typesAList)

cattedCoors=[];

for j = 1: length(MOLECULES)
    cattedCoors = cat(1,cattedCoors,MOLECULES(j).MOLCOORS);
end
COORDINATES = cattedCoors/LATTICE;
saveded = COORDINATES;
newCoords = zeros(0,3);
for m = 1: length(atomType)
  s = find(typesAList == atomType(m));
  newCoords = cat(1,newCoords, COORDINATES(s,:));
end
COORDINATES = newCoords;

