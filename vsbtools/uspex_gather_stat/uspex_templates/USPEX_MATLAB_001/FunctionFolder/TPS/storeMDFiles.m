function storeMDFiles(code, whichInd)

global ORG_STRUC
global TPS_STRUC

%
% We are now at the Calculation Folder
%
MDRestartFile = [' ' ORG_STRUC.MDrestartFile ' '];
aimFile       = [' ' ORG_STRUC.restartFile '_current' TPS_STRUC.POPULATION(whichInd).aim ' > convertRestart.log '];

if      code==4 % lammps
    %-- For Binay Lammps restart file, obsolete
    %
    %     if exist([ORG_STRUC.USPEXPath,'/FunctionFolder/TPS/restart2data'], 'file')
    %         convetor      = [ORG_STRUC.USPEXPath,'/FunctionFolder/TPS/./restart2data  ' ];
    %         [status, msg] = unix([convetor, MDRestartFile, aimFile]);
    %         if status > 0
    %             disp('Prolem met in converting the lammps binary restart file, please check ...');
    %             error(msg);
    %         end
    %     elseif exist(['./restart2data'], 'file')
    %         convetor      = ['./restart2data  ' ];
    %         [status, msg] = unix([convetor, MDRestartFile, aimFile]);
    %         if status > 0
    %             disp('Prolem met in converting the lammps binary restart file, please check ...');
    %             disp('Please visit: http://lammps.sandia.gov/doc/Section_tools.html#restart');
    %             error(msg);
    %         end
    %     else
    %         error('Cannot find restart2data for coverting the lammps binary restart file...');
    %     end
    [status, msg] = unix(['cp ' MDRestartFile ' ' aimFile]);
    if status > 0
        disp('Prolem met in copying the lammps restart file, please check ...');
        error(msg);
    end
    
elseif code==7 % CP2K
    
    
end

%  the restart file
restartFile = [ORG_STRUC.restartFile '_current' TPS_STRUC.POPULATION(whichInd).aim];
aimDiretory = [ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/'];
copyMDRestartFile(restartFile, aimDiretory);

% the trajectory file
restartFile = ORG_STRUC.trajectoryFile;
aimDiretory = [ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/',...
    ORG_STRUC.trajectoryFile '_current' TPS_STRUC.POPULATION(whichInd).aim];
copyMDRestartFile(restartFile, aimDiretory);

