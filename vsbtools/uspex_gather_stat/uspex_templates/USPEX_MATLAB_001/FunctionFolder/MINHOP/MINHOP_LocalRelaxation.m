function MINHOP_LocalRelaxation()

global ORG_STRUC
global POP_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     THE ACTUAL ALGORITHM    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while 1 % this eternal cycle is needed for non parallel calculations, parallelized one will /break out of it
         MINHOP_ReadJobs();
         MINHOP_SubmitJobs();
         if ORG_STRUC.platform > 0
           if sum([POP_STRUC.POPULATION(:).Done]) ~= length(POP_STRUC.POPULATION)
             delete('still_running');
             fclose ('all');
             quit
           else
             break; % break out of eternal cycle
           end
         else
           if sum([POP_STRUC.POPULATION(:).Done]) == length(POP_STRUC.POPULATION)
             break; % break out of eternal cycle for non parallel calculations 
           end
         end
    end   % end of 'eternal' cycle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  END MAIN LOOPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
