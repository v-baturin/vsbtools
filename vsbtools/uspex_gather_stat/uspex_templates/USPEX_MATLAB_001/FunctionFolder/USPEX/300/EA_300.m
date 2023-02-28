function EA_300()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global ANTISEEDS
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

ToStop = 0;

while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5
    %if ORG_STRUC.fixRndSeed > 0
    %	rng( ORG_STRUC.fixRndSeed+POP_STRUC.generation, 'twister' );
    %end
   
    if ORG_STRUC.currentGenDone == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LocalRelaxation();
        disp('status = Local optimization finished')
        disp(' ')
        fprintf(fp, [alignLine('Local optimization finished') '\n']);
        fprintf(fp, '\n');
        fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', ...
                                    POP_STRUC.generation) ) '\n'] );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 3:  Fitness                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SuccessRate();  %count how many structures were sucessfully done
        if ~((ORG_STRUC.pickedUP == 1) && ...
           (ORG_STRUC.pickUpGen == POP_STRUC.generation))

            fitness = CalcFitness_300();
            Correlation(fitness);
            fprintf(fp, [alignLine(sprintf(...
            'Correlation coefficient = %.4f',...
            POP_STRUC.correlation_coefficient) ) '\n'] );

            fitness = AntiSeedsCorrection(fitness);
            if ORG_STRUC.paretoRanking ~= 0
                for i=1:length(POP_STRUC.POPULATION)
                    Energy(i) = POP_STRUC.POPULATION(i).Enthalpies(end)...
                                /sum(POP_STRUC.POPULATION(i).numIons);
                end
		MINENERGY = min(Energy);
		Energy = Energy - MINENERGY;
                fitness = ParetoFitnessRankingCorrection(fitness, Energy);
            else
                fitness = FitnessRankingCorrection(fitness);
            end
            compositionStatistic(ORG_STRUC.firstGeneSplitAll);
            WriteGenerationOutput_300(fitness);
            if(ORG_STRUC.mlPeerFilter > 0 && ...
               ORG_STRUC.mlMinPopulation <= length(USPEX_STRUC.POPULATION))
                mlFilterTrain();
            end
        end

        ORG_STRUC.currentGenDone = 1;  %to avoid we calculate fitness again
        safesave([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC);
        safesave([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ToStop = 0;
    if ~((ORG_STRUC.pickedUP == 1) && ...
      (ORG_STRUC.pickUpGen == POP_STRUC.generation))
      [ToStop, Message] = StopRun(POP_STRUC.fitness, POP_STRUC.generation,...
      ORG_STRUC.numGenerations, ORG_STRUC.stopCrit, ORG_STRUC.stopFitness,...
      ORG_STRUC.fitLimit, ORG_STRUC.opt_sign * ORG_STRUC.optType);
    end
    if ToStop == 1
        fprintf(fp,'%30s \n', Message);
        Finish();
    else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        waitForCOPEX();

        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp,'\n');

        WriteGenerationBackup();
        fraction_content = update_STUFF('INPUT.txt');
        fprintf(fp, [alignLine('VARIATION OPERATORS') '\n']);
        fprintf(fp, '    Keep %2d percent to produce next generation\n',...
                                             ORG_STRUC.bestFrac*100);
        fprintf(fp, fraction_content);
        fprintf(fp,'\n');

        %determine appropriate coordinating numbers for topological random
        coordinating_number_analisis();
        CreateCalcFolder();

        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %make same struct
        OFF_STRUC.POPULATION(1) = QuickStart(OFF_STRUC.POPULATION);

        % A local flag for performing mlPeerFilter
        do_ml_filtering = false;

        if isfield(ORG_STRUC, 'mlPeerFilter') && (ORG_STRUC.mlPeerFilter > 0)
            min_generation = ceil(ORG_STRUC.mlMinPopulation /...
                                  ORG_STRUC.populationSize);
            if (POP_STRUC.generation >= min_generation)
                do_ml_filtering = true;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Operation = {'Heredity_300', 'Random_300', 'RandTop_300', ...
        'Permutation_300', 'LatticeMutation_300', 'SoftModeMutation_300'};
        Num_Opera = [ORG_STRUC.fracGene, ORG_STRUC.fracRand, ...
                  ORG_STRUC.fracRandTop, ORG_STRUC.fracPerm, ...
                  ORG_STRUC.fracLatMut,  ORG_STRUC.fracAtomsMut];
        Num_Opera = round(ORG_STRUC.populationSize * Num_Opera);
        count = 1;
        count_end = ORG_STRUC.populationSize;
        if do_ml_filtering
            count_end = ORG_STRUC.populationSize * ORG_STRUC.mlFilterRatio;
            Num_Opera = Num_Opera * ORG_STRUC.mlFilterRatio;
        end

        for i = 1 : length(Num_Opera)
            for j = 1 : Num_Opera(i)
                if count <= count_end
                    eval([Operation{i} '(' num2str(count) ')'])
                    %SOFTMutation might give two structures);
                    count = length(OFF_STRUC.POPULATION) + 1; 
                end
            end
        end
        if exist('spinOperation')
            spinOperation(count-1);   %NEED TO CHECK.........
        end

        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme');

        if do_ml_filtering
            OFF_STRUC.POPULATION = mlFilterPredict(OFF_STRUC.POPULATION,...
                                                 ORG_STRUC.populationSize);
        end

        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC.generation = POP_STRUC.generation + 1;
        fprintf(fp, [alignLine('Variation operators applied') '\n']);

        table = [strread(num2str(Num_Opera),'%s')'; Operation];
        fprintf(fp,' %10s structures produced by %s \n', table{:});

        addon_diff = KeepBestStructures();
        fprintf(fp,' %10d structures kept from previous gen\n',addon_diff);
        [addon_copex,OFF_STRUC] = importCOPEXStructures...
                                  (OFF_STRUC, POP_STRUC.generation);
        fprintf(fp,' %10d structures from other USPEX Calcs\n', addon_copex);

        POP_STRUC = OFF_STRUC;
%       num = pick_Seeds();
        num = pickUpSeeds();
        fprintf(fp,' %10d structures added from Seeds/POSCARS\n', num);
        fprintf(fp, [alignLine('Proceeding to new generation ') '\n']);
        fprintf(fp,' %10d parallel calculations running simutaneously\n',...
                                                 ORG_STRUC.numParallelCalcs);
        fprintf(fp, '\n');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fp, [alignLine( sprintf('Generation %d', POP_STRUC.generation) ) '\n'] );
        fprintf(fp,'  ID   Origin      Composition  Enthalpy(eV)');
        fprintf(fp,'  Volume(A^3)  KPOINTS  SYMMETRY\n');

        WriteGenerationStart();
        Start_POP();

        ORG_STRUC.currentGenDone = 0; %important to avoid fitness calculation
        if ( ORG_STRUC.pluginType > 0 ) && ...
           ( mod( POP_STRUC.generation, ORG_STRUC.pluginType) == 0 )
            ORG_STRUC.startNextGen = 0;
        else
            ORG_STRUC.startNextGen = 1;
        end
        
        safesave('Current_POP.mat', POP_STRUC);
        safesave('Current_ORG.mat', ORG_STRUC);
        disp(' ');
        disp('New generation produced');
        
    end
end
fclose(fp);
