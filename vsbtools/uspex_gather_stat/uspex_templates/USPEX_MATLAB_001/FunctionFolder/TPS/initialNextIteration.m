function initialNextIteration()

global ORG_STRUC
global TPS_STRUC


OFF_STRUC = createTPSStruct(ORG_STRUC.popSize);

disp(' ')
disp('<-> Initial next TPS iteration ...')

accept   = TPS_STRUC.shifter.accept;
direction= TPS_STRUC.direction;
successCounter = TPS_STRUC.successCounter;
maxFail  = ORG_STRUC.maxErrors;
%
% update the restartfile for next iteration
%
howCome = 'Shooter';
if TPS_STRUC.success == 1
    if TPS_STRUC.shifter.accept
        %=====================
        disp(' ')
        disp(' <-> Backup the last restart file')
        restartFile = ORG_STRUC.restartFile;
        aimDiretory = [ORG_STRUC.restartFile,'_last',  TPS_STRUC.direction(3)];
        copyMDRestartFile(restartFile, aimDiretory);
        
        %=====================
        disp(' ')
        disp(' <-> Update the new restart file');
        restartFile = [ORG_STRUC.restartFile,'_current',TPS_STRUC.direction(3)];
        aimDiretory = ORG_STRUC.restartFile;
        copyMDRestartFile(restartFile, aimDiretory);
        howCome = 'Shifter';

    else
        disp(' ')
        disp(' <-> NOT update the restart file');
    end
else
    %======================================================================================%
    %                                                                                      %
    %      ALL below lines are just for testing when droping in a metastable phase         %
    %                                                                                      %
    %======================================================================================%
    if     TPS_STRUC.successCounter.TPS <= -maxFail*2 && mod(TPS_STRUC.successCounter.TPS, maxFail) == 0;
        disp(' ')
        disp(' <-> Recovery from last opppsite direction restart file');
        disp('   -> Initial the shooter amplitude ...')
        restartFile = [ORG_STRUC.restartFile,'_last',TPS_STRUC.direction(1)];
        aimDiretory = ORG_STRUC.restartFile;
        copyMDRestartFile(restartFile, aimDiretory);
        howCome = 'Revert2';
        TPS_STRUC.amplitudeA2B = ORG_STRUC.amplitudeShoot(1);
        TPS_STRUC.amplitudeB2A = ORG_STRUC.amplitudeShoot(2);

    elseif TPS_STRUC.successCounter.TPS == -maxFail 
        disp(' ')
        disp(' <-> Recovery from last restart file');
        restartFile = [ORG_STRUC.restartFile,'_last',TPS_STRUC.direction(3)];
        aimDiretory = ORG_STRUC.restartFile;
        copyMDRestartFile(restartFile, aimDiretory);
        howCome = 'Revert';
        if TPS_STRUC.direction(3)=='A'
            if TPS_STRUC.amplitudeB2A < ORG_STRUC.amplitudeShoot(2)
                TPS_STRUC.amplitudeB2A = TPS_STRUC.amplitudeB2A*(TPS_STRUC.magnitudeB2A^(maxFail*1.5));
            else
                TPS_STRUC.amplitudeB2A = ORG_STRUC.amplitudeShoot(2);
            end
        else
            if TPS_STRUC.amplitudeA2B < ORG_STRUC.amplitudeShoot(1)
                TPS_STRUC.amplitudeA2B = TPS_STRUC.amplitudeA2B*(TPS_STRUC.magnitudeA2B^(maxFail*1.5));
            else
                TPS_STRUC.amplitudeA2B = ORG_STRUC.amplitudeShoot(1);
            end
        end
        
    else
        disp('The MD trajectory meets some problem. Please check...')
    end
    %======================================================================================%
    %                                                                                      %
    %                                          End 
    %                                                                                      %
    %======================================================================================%
end


%
% Initial the TPS_STRUC
%
if accept == 1
    disp([' -> Current TPS trajectory(' direction ') succeeded'])
    OFF_STRUC = initialOFFSTRUC(OFF_STRUC, TPS_STRUC, accept, successCounter, howCome);
else
    disp([' -> Current TPS trajectory(' direction ') failed']);
    OFF_STRUC = initialOFFSTRUC(OFF_STRUC, TPS_STRUC, accept, successCounter, howCome);
end


% new iteration
TPS_STRUC = OFF_STRUC;


cd([ORG_STRUC.homePath, '/', ORG_STRUC.resFolder ]);
makeIterationFolder();



disp('<-> Initializaiotn DONE.')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OFF_STRUC = initialOFFSTRUC(OFF_STRUC, TPS_STRUC, accept, successCounter, howCome)


OFF_STRUC.iteration = TPS_STRUC.iteration+1;
disp([' -> Next Iiteration Num : ' num2str(OFF_STRUC.iteration)] );


% Define the opposite direction
if strcmp(TPS_STRUC.direction, 'A2B')
    oppoDiret = 'B2A';
else
    oppoDiret = 'A2B';
end

if accept
    OFF_STRUC.direction = oppoDiret;
else
    OFF_STRUC.direction = TPS_STRUC.direction;
end

OFF_STRUC.success   = 0;
OFF_STRUC.howCome   = howCome;

OFF_STRUC.AHT = [];
OFF_STRUC.BHT = [];
OFF_STRUC.lastAHT = TPS_STRUC.AHT;
OFF_STRUC.lastBHT = TPS_STRUC.BHT;
OFF_STRUC.successCounter = TPS_STRUC.successCounter;

OFF_STRUC.magnitudeA2B = TPS_STRUC.magnitudeA2B;
OFF_STRUC.magnitudeB2A = TPS_STRUC.magnitudeB2A;
OFF_STRUC.amplitudeA2B = TPS_STRUC.amplitudeA2B;
OFF_STRUC.amplitudeB2A = TPS_STRUC.amplitudeB2A;
if accept
    disp(' -> Current TPS trajectory accepted')
    disp([' -> Next direction : ' OFF_STRUC.direction ]);
    disp([' -> Determine shooter amplitude for ', TPS_STRUC.direction]);
    if OFF_STRUC.direction(end) == 'A'
        if successCounter.B2A > 5
            OFF_STRUC.amplitudeB2A = OFF_STRUC.amplitudeB2A*OFF_STRUC.magnitudeB2A;
        end
    else
        if successCounter.A2B > 5
            OFF_STRUC.amplitudeA2B = OFF_STRUC.amplitudeA2B*OFF_STRUC.magnitudeA2B;
        end
    end
else
    disp(' -> Current TPS trajectory rejected')
    disp([' -> Next direction : ' OFF_STRUC.direction ]);
    disp([' -> Decrease shooter amplitude for ', TPS_STRUC.direction]);
    if OFF_STRUC.direction(end) == 'A'
        OFF_STRUC.amplitudeB2A = OFF_STRUC.amplitudeB2A/OFF_STRUC.magnitudeB2A;
    else
        OFF_STRUC.amplitudeA2B = OFF_STRUC.amplitudeA2B/OFF_STRUC.magnitudeA2B;
    end
end


OFF_STRUC.POPULATION(1).aim = OFF_STRUC.direction(1);
OFF_STRUC.POPULATION(2).aim = OFF_STRUC.direction(3);

disp(['  -> CalcFold1 Aim : ' OFF_STRUC.POPULATION(1).aim] );
if OFF_STRUC.POPULATION(1).aim == 'A'
    disp(['  -> Amplitude B2A : ' num2str(OFF_STRUC.amplitudeB2A)] );
else
    disp(['  -> Amplitude A2B : ' num2str(OFF_STRUC.amplitudeA2B)] );
end
disp(['  -> CalcFold2 Aim : ' OFF_STRUC.POPULATION(2).aim] );
if OFF_STRUC.POPULATION(2).aim == 'B'
    disp(['  -> Amplitude A2B : ' num2str(OFF_STRUC.amplitudeA2B)] );
else
    disp(['  -> Amplitude B2A : ' num2str(OFF_STRUC.amplitudeB2A)] );
end

OFF_STRUC.lastrun.success   = TPS_STRUC.success;
OFF_STRUC.lastrun.direction = TPS_STRUC.direction;



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
