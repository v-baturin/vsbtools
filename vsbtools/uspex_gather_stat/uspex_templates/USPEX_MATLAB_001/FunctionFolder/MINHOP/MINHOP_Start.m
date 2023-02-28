function MINHOP_Start()

%Add USPEX_STRUC
%Lastly updated by Qiang Zhu, 2014/02/18
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

% initialize the random number generator
initialRandomGenerator();


disp('  ');
while 1
    stillRunningInfo();
    [homePath, USPEXPath] = workingPath();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% To get POP and ORG_STRUC, etc......    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if exist('Current_POP.mat')
        load Current_POP.mat
        load Current_ORG.mat
        %%%%%%%%%
        verifyWorkingPath(homePath, USPEXPath);
	%%%%%%%%%
        cd (POP_STRUC.resFolder)
        load USPEX.mat
        cd ..
        
    else
        inputFile = 'INPUT.txt';
        MINHOP_CreateORGStruc(inputFile);
        unixCmd (['cp ' inputFile ' ' ORG_STRUC.resFolder '/Parameters.txt']);
        unixCmd (['cp POSCAR_1 '      ORG_STRUC.resFolder '/']);

        %ORG_STRUC.pickedUP = 0;
        CreateCalcFolder(1);
        if ORG_STRUC.pickUpYN
            PickUp();
        else
            MINHOP_Initialize();
        end
        CreateCalcFolder();
        
    end
    
    safesave ('Current_POP.mat', POP_STRUC)
    safesave ('Current_ORG.mat', ORG_STRUC)

    MINHOP();
    
    if ORG_STRUC.repeatForStatistics > 1
        RepeatRun();
    else
        USPEXDone();
    end
end
