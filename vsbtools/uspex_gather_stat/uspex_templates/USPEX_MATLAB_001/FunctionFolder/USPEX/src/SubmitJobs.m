function SubmitJobs(ID_next, indic)
global ORG_STRUC
global POP_STRUC

DO_NOW = ID_next;
N_Step = length(ORG_STRUC.abinitioCode);

POP_STRUC.POPULATION(DO_NOW).Folder = indic;
POP_STRUC.POPULATION(DO_NOW).ToDo   = 0;
Step        = POP_STRUC.POPULATION(DO_NOW).Step;
LATTICE     = POP_STRUC.POPULATION(DO_NOW).LATTICE;
COORDINATES = POP_STRUC.POPULATION(DO_NOW).COORDINATES;

if (ORG_STRUC.dimension == 3) && (Step > 1) && (Step < N_Step - 1)
    LATTICE     = perturbCell(LATTICE);
    COORDINATES = perturbCoords(COORDINATES, LATTICE,1); 
    POP_STRUC.POPULATION(DO_NOW).LATTICE     = LATTICE;
    POP_STRUC.POPULATION(DO_NOW).COORDINATES = COORDINATES;
end

if  Step > N_Step % structures from Best
    POP_STRUC.POPULATION(DO_NOW).JobID = 0.01;
elseif ORG_STRUC.abinitioCode(Step) == 0   % no optimization at all! (used in order optimization)
    POP_STRUC.POPULATION(DO_NOW).JobID = 0.02;
else
    if ORG_STRUC.dimension == 0
       lat = POP_STRUC.POPULATION(DO_NOW).LATTICE;
       candidate = POP_STRUC.POPULATION(DO_NOW).COORDINATES;
       [lat, candidate] = reduce_Cluster(lat, candidate);
       [lat, candidate] = makeCluster(lat, candidate, ORG_STRUC.vacuumSize(Step));
       [candidate] = moveCluster(lat, candidate);
       POP_STRUC.POPULATION(DO_NOW).LATTICE = lat;
       POP_STRUC.POPULATION(DO_NOW).COORDINATES = candidate;
    end
    cd (['CalcFold' num2str(indic)])
    Write_AbinitCode(ORG_STRUC.abinitioCode(Step), DO_NOW);
    POP_STRUC.POPULATION(DO_NOW).JobID = submitJob(DO_NOW);
    cd ..
end

safesave ('Current_POP.mat', POP_STRUC)
