function TPS()

global ORG_STRUC
global TPS_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%                     TPS   ALGORITHM             %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (TPS_STRUC.iteration < ORG_STRUC.numIterations + 0.5)
    cd (ORG_STRUC.homePath);
    
    while 1 % this eternal cycle is needed for non parallel calculations, parallelized one will break out of it
        
        %KillJobs_TPS();
        ReadJobs_TPS();
        SubmitJobs_TPS();
        
        if (ORG_STRUC.platform > 0) || (ORG_STRUC.numParallelCalcs > 1)
            if sum([TPS_STRUC.POPULATION(:).Done])~= length(TPS_STRUC.POPULATION)
                quitUSPEX();
            else
                break;
            end
        else
            if sum([TPS_STRUC.POPULATION(:).Done]) == length(TPS_STRUC.POPULATION)
                break; % break out of eternal cycle when non parallel calculations finished for given generation
            end
        end
    end   % end of 'eternal' cycle
    
    
    checkTPSSuccess();
    runShifter();
    
    saveTPSFiles();
    WriteTPSOutput();
    plotTPSFiguers();
    %
    initialNextIteration();
    runShooter();

    
    cd (ORG_STRUC.homePath);
    safesave ('Current_ORG.mat', ORG_STRUC)
    safesave ('Current_TPS.mat', TPS_STRUC)
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function quitUSPEX()

unixCmd('rm still_running');
quit();



