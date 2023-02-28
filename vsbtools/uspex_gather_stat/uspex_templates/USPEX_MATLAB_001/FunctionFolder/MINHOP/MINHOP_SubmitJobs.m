function MINHOP_SubmitJobs()
global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    DO_NOW = 0;
    if isempty (whichInd)
        stillLeft = find([POP_STRUC.POPULATION(:).ToDo]);
        if ~isempty(stillLeft)
           DO_NOW = stillLeft(1);
           POP_STRUC.POPULATION(DO_NOW).Folder = indic;
           POP_STRUC.POPULATION(DO_NOW).ToDo = 0;
        end
    elseif POP_STRUC.POPULATION(whichInd).JobID == 0
        DO_NOW = whichInd;
    end
    if DO_NOW
        Step = POP_STRUC.POPULATION(DO_NOW).Step;
        cd (['CalcFold' num2str(indic)])
        Write_AbinitCode(ORG_STRUC.abinitioCode(Step), DO_NOW);
        POP_STRUC.POPULATION(DO_NOW).JobID = submitJob(DO_NOW);
        cd ..
        safesave ('Current_POP.mat', POP_STRUC)
    end
end
