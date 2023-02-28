function Start_PSO()

global ORG_STRUC
global POP_STRUC
global PSO_STRUC
global ANTISEEDS

% initialize the random number generator
initialRandomGenerator();


disp('  ');

while 1
    stillRunningInfo();
    [homePath, USPEXPath] = workingPath();

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% To get POP and ORG_STRUC, etc......    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if exist('Current_POP.mat')
        load Current_POP.mat;
        load Current_ORG.mat;
        verifyWorkingPath(homePath, USPEXPath);
        %%%%%%%%%
        cd (POP_STRUC.resFolder)
        load PSO.mat

        cd ..

        if exist('ANTISEEDS.mat')
            load ANTISEEDS.mat;
        else
            pickAntiSeeds();
        end

    else
        inputFile = 'INPUT.txt';
        createORGStruc(inputFile);
        unixCmd (['cp ' inputFile ' ' ORG_STRUC.resFolder '/Parameters.txt']);
        %ORG_STRUC.pickedUP = 0;

        % create CalcFoldTemp Folder
        CreateCalcFolder(1);

        if ORG_STRUC.pickUpYN
            PickUp();
        else
            Initialize_PSO();
        end
        CreateCalcFolder();
    end

    %    safesave('Current_POP.mat', 'POP_STRUC')
    %    safesave('Current_ORG.mat', 'ORG_STRUC')
    save('Current_POP.mat', 'POP_STRUC');
    save('Current_ORG.mat', 'ORG_STRUC');

    path(path,[homePath '/FunctionFolder/PSO']);
    PSO();

    if ORG_STRUC.repeatForStatistics > 1
        RepeatRun();
    else
        USPEXDone();
    end
end
