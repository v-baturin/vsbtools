function Start()

%Add USPEX_STRUC
%Lastly updated by Qiang Zhu, 2014/02/18
global ORG_STRUC
global POP_STRUC
global POOL_STRUC
global ANTISEEDS
global USPEX_STRUC
global FORCE_STRUC

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

     	checkCompatibility();
        verifyWorkingPath(homePath, USPEXPath);

		%%%%%%%%%
        cd (POP_STRUC.resFolder)
        load USPEX.mat
        if exist('FORCE.mat')
           load FORCE.mat
        end
        if exist('POOL.mat')
            load POOL.mat
        end
        cd ..
        
        if exist('ANTISEEDS.mat')
            load ANTISEEDS.mat
        else
            pickAntiSeeds();
        end
		
		if ORG_STRUC.fixRndSeed > 0 
			rng( ORG_STRUC.fixRndSeed+POP_STRUC.generation, 'twister' );
		end
    else
        inputFile = 'INPUT.txt';
        createORGStruc(inputFile);

        unixCmd(['cp ' inputFile ' '    ORG_STRUC.resFolder '/Parameters.txt']);
        unixCmd(['cp MOL_* '            ORG_STRUC.resFolder '/']);
        unixCmd(['cp POSCAR_SUBSTRATE ' ORG_STRUC.resFolder '/']);

        % Preprocess data before EA selection:
        %if ((ORG_STRUC.dimension==3) && (ORG_STRUC.molecule==1)) || ...
        %   ((ORG_STRUC.dimension==1) && (ORG_STRUC.molecule==1) && (ORG_STRUC.varcomp==0))
        %   [a,b]=unix (['cp MOL_* ' ORG_STRUC.resFolder '/']);
        %elseif (ORG_STRUC.dimension==2) && (ORG_STRUC.molecule==0)
        %   [a,b]=unix (['cp POSCAR_SUBSTRATE ' Current_ORG_STRUC.resFolder '/']);
        %end


        path(path,[USPEXPath '/FunctionFolder/USPEX/' calcTypeStr]);
        eval(['PreCheck_' calcTypeStr '()']);

        %ORG_STRUC.pickedUP = 0;
        
        % create CalcFoldTemp Folder
        CreateCalcFolder(1);
        if ORG_STRUC.pickUpYN
            PickUp();
        else
            Initialize();
        end
        CreateCalcFolder();
    end
   
    safesave ('Current_POP.mat', POP_STRUC)
    safesave ('Current_ORG.mat', ORG_STRUC)
    

    % EA selection part:

    path(path,[USPEXPath '/FunctionFolder/USPEX/' calcTypeStr]);
    eval(['EA_' calcTypeStr '()']);
    
        
    if ~isempty(ORG_STRUC.stopFitness)
        disp('Run statistics analysis:');
        statistics(ORG_STRUC.stopFitness);
    end

    if ORG_STRUC.repeatForStatistics > 1
        RepeatRun();
    else
        USPEXDone();
    end
    
end


