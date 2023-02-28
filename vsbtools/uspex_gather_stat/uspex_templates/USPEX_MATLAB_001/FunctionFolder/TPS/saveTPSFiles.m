function saveTPSFiles()


global ORG_STRUC
global TPS_STRUC


cd([ORG_STRUC.homePath, '/', ORG_STRUC.resFolder ]);
makeIterationFolder();

%--------------------------%
%  Update MD restart file  %
%--------------------------%
    
% store the restart files to iteration folder
aimDiretory = ['Iterations/iter-' num2str(TPS_STRUC.iteration),'/' ];

restartFile = ORG_STRUC.restartFile;
copyMDRestartFile(restartFile, aimDiretory);

restartFile = [ORG_STRUC.restartFile,'_currentA' ];
copyMDRestartFile(restartFile, aimDiretory);

restartFile = [ORG_STRUC.restartFile,'_currentB' ];
copyMDRestartFile(restartFile, aimDiretory);


%--------------------------%
%   Save shoot log files   %
%--------------------------%
disp(' ')
disp(' <-> Save shootor log files');
tempFolder = [ORG_STRUC.homePath,'/CalcFoldTemp/'];
aimDiretory = ['Iterations/iter-' num2str(TPS_STRUC.iteration),'/' ];
unixCmd(['cp ' tempFolder '/shoot.log   ' aimDiretory]);
unixCmd(['cp ' tempFolder '/input_shoot ' aimDiretory]);
unixCmd(['cp ' tempFolder '/lastrun     ' aimDiretory]);


%--------------------------%
% Save MD trajectory files %
%--------------------------%
disp(' ')
disp(' <-> Save the whole trajectory file');
resFolder   = [ORG_STRUC.homePath,'/',ORG_STRUC.resFolder];
aimFile = ['Iterations/iter-' num2str(TPS_STRUC.iteration),'/',ORG_STRUC.trajectoryFile ];
invsersPy   = [ORG_STRUC.USPEXPath,'/FunctionFolder/TPS/inverseXYZfile.py'];

sourceFile1 = [resFolder, '/', ORG_STRUC.trajectoryFile, '_current', TPS_STRUC.direction(1)];
sourceFile2 = [resFolder, '/', ORG_STRUC.trajectoryFile, '_current', TPS_STRUC.direction(3)];

python_uspex( invsersPy, [' ' sourceFile1 ' > ' aimFile ]);
[status, msg] = unix([ 'cat ' sourceFile2 ' >> ' aimFile]);
if status > 0
    disp(msg);
    error('Failed to creat the completet MD trajectory file.');
end


%--------------------------%
%   Save MSD & CN  files   %
%--------------------------%
calcFolder1 = [ORG_STRUC.homePath,'/CalcFold1/'];
calcFolder2 = [ORG_STRUC.homePath,'/CalcFold2/'];

try
        disp(' <-> Save MSD file ')
        unixCmd(['cp ' calcFolder1 '/*.msd    ' aimDiretory]);
        unixCmd(['cp ' calcFolder2 '/*.msd    ' aimDiretory]);
catch
        disp(['---Failed to copy ...'])
end

try
        disp(' <-> Save CRN file ')
        unixCmd(['cp ' calcFolder1 '/*.crn    ' aimDiretory]);
        unixCmd(['cp ' calcFolder2 '/*.crn    ' aimDiretory]);
catch
        disp(['---Failed to copy ...'])
end


% save the ORG_STRUC and TPS_STRUC file
cd(aimDiretory);
safesave ('ORG_STRUC.mat', ORG_STRUC);
safesave ('TPS_STRUC.mat', TPS_STRUC);
cd(ORG_STRUC.homePath);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
%
%
function makeIterationFolder()

global TPS_STRUC

if ~exist('Iterations', 'dir')
    mkdir('Iterations');
end

cd('Iterations');
if ~exist(['iter-' num2str(TPS_STRUC.iteration) ], 'dir')
    mkdir(['iter-' num2str(TPS_STRUC.iteration) ]);
end
cd ..
