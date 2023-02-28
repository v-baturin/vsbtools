function SubmitJobs_TPS()


global ORG_STRUC
global TPS_STRUC
global POP_STRUC  % to have submitJobs compitiable.

disp(' ');
disp('=== Submit jobs here :')
for indic = 1:ORG_STRUC.numParallelCalcs
    whichInd = find([TPS_STRUC.POPULATION(:).Folder]==indic);
    DO_NOW = 0;
    if isempty (whichInd)
        stillLeft = find([TPS_STRUC.POPULATION(:).ToDo]);
        if ~isempty(stillLeft)
            DO_NOW = stillLeft(1);
            TPS_STRUC.POPULATION(DO_NOW).Folder = indic;
            TPS_STRUC.POPULATION(DO_NOW).ToDo = 0;
        end
    elseif TPS_STRUC.POPULATION(whichInd).JobID == 0
        DO_NOW = whichInd;
    end
    
    if DO_NOW
        Step = TPS_STRUC.POPULATION(DO_NOW).Step;
        if Step > length ([ORG_STRUC.abinitioCode]) % structures from Best
            TPS_STRUC.POPULATION(DO_NOW).JobID = 0.01;
        elseif ORG_STRUC.abinitioCode(Step) == 0   % no optimization at all! (used in order optimization)
            TPS_STRUC.POPULATION(DO_NOW).JobID = 0.02;
        else
            cd ([ORG_STRUC.homePath '/CalcFold' num2str(indic)]);
            TPSAbinitWriting(ORG_STRUC.abinitioCode(Step), DO_NOW);
            POP_STRUC = TPS_STRUC;
            TPS_STRUC.POPULATION(DO_NOW).JobID = submitJob(DO_NOW);
                        
            cd (ORG_STRUC.homePath)
        end
        safesave ('Current_TPS.mat', TPS_STRUC)
    end
end

cd (ORG_STRUC.homePath)

disp('=== Job submission DONE.')

