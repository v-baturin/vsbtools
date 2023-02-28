function NEB_Start()

% to start a VCNEB code

global ORG_STRUC
global POP_STRUC

disp('  ');
while 1
    
    checkUSPEXDone();
    stillRunningInfo();
    
    [homePath, USPEXPath] = workingPath();
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% To get POP and ORG_STRUC, etc......    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('Current_POP.mat')
        load Current_POP.mat
        load Current_ORG.mat
        
        %checkCompatibility();
        %%%%%%%%%
        verifyWorkingPath(homePath, USPEXPath);
        %%%%%%%%%
        cd (POP_STRUC.resFolder)
        %load USPEX.mat
    else
        inputFile = 'INPUT.txt';
        
        createORGStruc(inputFile);
        
        unixCmd (['cp ' inputFile ' ' ORG_STRUC.resFolder '/Parameters.txt']);
        unixCmd (['cp POSCAR_1 '      ORG_STRUC.resFolder '/initialImages']);
        %ORG_STRUC.pickedUP = 0;
        
        % create CalcFoldTemp Folder
        CreateCalcFolder(1);
        if ORG_STRUC.pickUpYN
            NEB_PickUp();
        else
            NEB_initialize_POP_STRUC();
        end
        % create CalcFoldID Folders
        CreateCalcFolder();
        POP_STRUC.resFolder = ORG_STRUC.resFolder;
    end
    
    safesave ([ORG_STRUC.homePath, '/Current_POP.mat'], POP_STRUC)
    safesave ([ORG_STRUC.homePath, '/Current_ORG.mat'], ORG_STRUC)
    
    % if ORG_STRUC.
    VCNEB();
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

