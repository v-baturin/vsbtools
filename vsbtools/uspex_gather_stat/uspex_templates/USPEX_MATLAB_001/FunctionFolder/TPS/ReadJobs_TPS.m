function ReadJobs_TPS()

global ORG_STRUC
global TPS_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    cd (ORG_STRUC.homePath)
    
    % check whether calculations are still performed and read out the output if they are done
    whichInd = find([TPS_STRUC.POPULATION(:).Folder]==indic);
    if ~isempty (whichInd)
        disp(' ');
        disp(['<-> MD calculation in CalcFold' num2str(whichInd) ' .......']);
        Step = TPS_STRUC.POPULATION(whichInd).Step;
        if TPS_STRUC.POPULATION(whichInd).JobID
            if (ORG_STRUC.platform > 0) || (ORG_STRUC.numParallelCalcs > 1)
                disp(['JobID=' num2str(TPS_STRUC.POPULATION(whichInd).JobID) ]);
            end
            POP_STRUC = TPS_STRUC;
            doneOr = checkStatusC(whichInd);
            cd ([ORG_STRUC.homePath '/CalcFold' num2str(indic)]);
            if doneOr
                Error = TPSAbinitReading(ORG_STRUC.abinitioCode(Step), whichInd);
                disp(' --> Job Done.')
                if Error == 0
                    disp('  -> No error, store the MD restart file.')
                    % TPSAbinitClean(ORG_STRUC.abinitioCode(Step));    % clean the files in CalcFold when done
                    storeMDFiles(ORG_STRUC.abinitioCode(Step), whichInd);
                    disp('  -> MD restart file stored.')
                    TPS_STRUC.POPULATION(whichInd).Done = 1;
                    TPS_STRUC.bodyCount = TPS_STRUC.bodyCount + 1;
                    TPS_STRUC.DoneOrder(whichInd) = TPS_STRUC.bodyCount;
                    TPS_STRUC.POPULATION(whichInd).Folder=0;
                else
                    disp('  !-> Error found, job will be ran again!')
                    TPS_STRUC.POPULATION(whichInd).ToDo=1;
                end
                cd (ORG_STRUC.homePath)
                safesave ('Current_TPS.mat', TPS_STRUC);
                TPS_STRUC.POPULATION(whichInd).JobID = 0;
            else
                disp(' --> Not finshed yet.')
            end            
        end
        disp(['<-> Job Checking Done']);
    end
    cd (ORG_STRUC.homePath);
    safesave ('Current_TPS.mat', TPS_STRUC);
    
end
cd (ORG_STRUC.homePath)
