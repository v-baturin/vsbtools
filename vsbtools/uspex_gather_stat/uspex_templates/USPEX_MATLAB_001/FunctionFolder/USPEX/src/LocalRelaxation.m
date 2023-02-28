function LocalRelaxation()

%This rountine is used to handle Jobs in local optmization
%1, check all pending jobs (ReadJobs**.m)
%2, submit new jobs (SubmitJobs_**.m) if needed

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

dimension  = ORG_STRUC.dimension;
molecule   = ORG_STRUC.molecule;
varcomp    = ORG_STRUC.varcomp;
N_Parallel = ORG_STRUC.numParallelCalcs;
platform   = ORG_STRUC.platform;
N_popsize  = length(POP_STRUC.POPULATION);

calcType   = [num2str(dimension) '' num2str(molecule) '' num2str(varcomp)];
if str2num(calcType) < 0
    tmp = -1*str2num(calcType);
    calcType = ['M' num2str(tmp)];
end

%The first eternal cycle is needed for non parallel calculations,
%parallelized one will /break out of it
while 1 
    for indic = 1:N_Parallel
        ID_now = POP_STRUC.CalcFold(indic);
        if (ID_now>0) && (ID_now<=N_popsize)
           Calc_max   = POP_STRUC.CalcFold_max;
           ID_next = ReadJobs(ID_now, Calc_max, indic);
           if (ID_next > 0) && (ID_next <= N_popsize)
              SubmitJobs(ID_next, indic);
           end
        end
    end

    if (platform > 0) || (N_Parallel > 1)
      if sum([POP_STRUC.POPULATION(:).Done])~= N_popsize
         delete('still_running');  %Remove still_running
         fclose ('all');           %exit Matlab
         quit
      else
         break; % break out of eternal cycle for parallel mode
      end
    else
       if sum([POP_STRUC.POPULATION(:).Done]) == N_popsize
           gen = POP_STRUC.generation;
           disp(['The end of generation ' num2str(gen)]);
           if ~exist('restart_matlab') && (1*fix(gen/1)==gen)
               safesave([ORG_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
               delete('still_running');
               unixCmd('echo "This file is printed when matlab should be reatarted"   > restart_matlab');
               unixCmd('echo "to prevent memory overflow and slowing down of matlab" >> restart_matlab');
               fclose ('all');           %exit Matlab
               quit               
           else
               delete('restart_matlab');
           end
           break; % break out of eternal cycle for non parallel mode
       end
    end
end

