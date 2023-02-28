function TPS_Start()

% to start a VCNEB code

global ORG_STRUC
global TPS_STRUC


disp('  ');
while 1
    stillRunningInfo();
    [homePath, USPEXPath] = workingPath();

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% To get POP and ORG_STRUC, etc......    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('Current_TPS.mat')
        load Current_TPS.mat
        load Current_ORG.mat
        
        %checkCompatibility();
        %%%%%%%%%
        verifyWorkingPath(homePath, USPEXPath);
        %%%%%%%%%
        cd (ORG_STRUC.resFolder)
        %load USPEX.mat
    else
        inputFile = 'INPUT.txt';
        
        createORGStruc(inputFile);
        unixCmd (['cp ' inputFile ' ' ORG_STRUC.resFolder '/Parameters.txt']);
        
        copyMDRestartFile([ORG_STRUC.homePath '/Seeds/' ORG_STRUC.restartFile],...
                          [ORG_STRUC.homePath '/'       ORG_STRUC.resFolder '/' ORG_STRUC.restartFile]);
        
        % create CalcFoldTemp Folder
        CreateCalcFolder(1);
        if ORG_STRUC.pickUpYN
            TPS_PickUp();
        else
            initialize_TPS_STRUC();
        end
        
        % create CalcFoldID Folders
        CreateCalcFolder();
        
        TPS_STRUC.resFolder = ORG_STRUC.resFolder;
    end
      
    safesave ([ORG_STRUC.homePath, '/Current_TPS.mat'], TPS_STRUC)
    safesave ([ORG_STRUC.homePath, '/Current_ORG.mat'], ORG_STRUC)
    
    % if ORG_STRUC.
    TPS();
    % else
    % tempertureRelxa();
    %end
    
    if ORG_STRUC.repeatForStatistics > 1
        RepeatRun();
    else
        USPEXDone();
    end
    
end


%%--------------------------------------
%%--------------------------------------

