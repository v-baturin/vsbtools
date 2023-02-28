function ReadJobs_NEB()

global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    cd (ORG_STRUC.homePath)
    
    % check whether calculations are still performed and read out the output if they are done
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    if ~isempty (whichInd)
        disp(' ');
        disp(['Image ' num2str(whichInd) '.......']);
        Step = POP_STRUC.POPULATION(whichInd).Step;
        if POP_STRUC.POPULATION(whichInd).JobID
            if (ORG_STRUC.platform > 0) || (ORG_STRUC.numParallelCalcs > 1)
                disp(['JobID=' num2str(POP_STRUC.POPULATION(whichInd).JobID) ]);
            end
            doneOr = checkStatusC(whichInd);
            if doneOr                
                cd ([ORG_STRUC.homePath '/CalcFold' num2str(indic)]);
                Error = Reading(ORG_STRUC.abinitioCode(Step),whichInd, indic);
                if POP_STRUC.POPULATION(whichInd).Error == 0
                    disp('-> Job Done Successfully');

                    POP_STRUC.POPULATION(whichInd).Done = 1;
                    POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
                    POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
                    POP_STRUC.POPULATION(whichInd).Folder=0;
                else
                    disp('-> Job Failed');

                    POP_STRUC.POPULATION(whichInd).ToDo=1;
                    POP_STRUC.POPULATION(whichInd).JobID = 0;
                    POP_STRUC.POPULATION(whichInd).Folder=0;
                end
                cd (ORG_STRUC.homePath)
                safesave ('Current_POP.mat', POP_STRUC);
                POP_STRUC.POPULATION(whichInd).JobID = 0;
            else
                 disp('-> Not finished yet');
            end
        end
    end
    cd (ORG_STRUC.homePath);
    safesave ('Current_POP.mat', POP_STRUC);
end
cd (ORG_STRUC.homePath)
