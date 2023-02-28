function initialize_TPS_STRUC()


global ORG_STRUC
global TPS_STRUC


%======== initial the trajectory information !  ==========%

ORG_STRUC.popSize = 2;


%
% ( A-> start point ->  B ) : Direction A2B or B
% ( B-> start point ->  A ) : Direction B2A or A
%
%
TPS_STRUC = createTPSStruct(ORG_STRUC.popSize);



TPS_STRUC.iteration = 0;
TPS_STRUC.direction = 'A2B';
TPS_STRUC.success   = 1;
TPS_STRUC.howCome   = 'inital';

TPS_STRUC.AHT = [];
TPS_STRUC.BHT = [];
TPS_STRUC.lastAHT = [1, 0, 0];
TPS_STRUC.lastBHT = [1, 0, 0];

TPS_STRUC.amplitudeA2B = ORG_STRUC.amplitudeShoot(1);
TPS_STRUC.magnitudeA2B = ORG_STRUC.magnitudeShoot(1);
TPS_STRUC.amplitudeB2A = ORG_STRUC.amplitudeShoot(2);
TPS_STRUC.magnitudeB2A = ORG_STRUC.magnitudeShoot(2);


TPS_STRUC.POPULATION(1).aim = TPS_STRUC.direction(1);
TPS_STRUC.POPULATION(2).aim = TPS_STRUC.direction(3);

TPS_STRUC.lastrun.success=1;
TPS_STRUC.lastrun.direction = 'B2A';

%======== Shoot the initial trajectory here !  ==========%

disp('')
disp('<-> Initial the first TPS trajectory ...')
%-- check the CalcFoldTemp, and access it
tmpFolder = [ORG_STRUC.homePath '/CalcFoldTemp' ];

if ~exist(tmpFolder, 'dir')
    mkdir(tmpFolder);
end
cd(tmpFolder);


%-- repeat the first trajectory to confirm it

delete(ORG_STRUC.restartFile); % delete the existing MD restart file.
restartFile = [ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/' ORG_STRUC.restartFile];
copyMDRestartFile(restartFile, 'ascii.new');

%runShooter();
%--


TPS_STRUC.success = 0;
TPS_STRUC.POPULATION(1).success = 0;
TPS_STRUC.POPULATION(2).success = 0;

disp('<-> Initialization DONE.')
