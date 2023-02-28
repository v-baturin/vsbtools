function EA_M300()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global USPEX_STRUC
global ANTISEEDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5

    if ORG_STRUC.currentGenDone == 0
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       LocalRelaxation();

       disp('status = Local optimization finished')
       disp(' ')
       fprintf(fp, [alignLine('-', 0) '\n']);
       fprintf(fp, [alignLine('Local optimization finished') '\n']);
       fprintf(fp, [alignLine('-', 0) '\n']);
       fprintf(fp, '\n');
       fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', ...
                                   POP_STRUC.generation) ) '\n'] );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 3:  Fitness          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SuccessRate();  %count how many structures were sucessfully done

        if ~((ORG_STRUC.pickedUP == 1) && ...
            (ORG_STRUC.pickUpGen == POP_STRUC.generation))
            fitness = CalcFitness_M300();
            cor_coefficient = Correlation(fitness);
            fprintf(fp, [alignLine(sprintf('Correlate coefficient = %.4f',...
                                               cor_coefficient) ) '\n'] );

            fitness = AntiSeedsCorrection(fitness);
            fitness = FitnessRankingCorrection(fitness);
            WriteGenerationOutput_M300(fitness);
        end
        ORG_STRUC.currentGenDone=1;  %to avoid calculate fitness again
        safesave ([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC)
        safesave ([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC)
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ToStop = 0;
    if ~((ORG_STRUC.pickedUP == 1) && ...
        (ORG_STRUC.pickUpGen == POP_STRUC.generation))
        [ToStop, Message] = StopRun(POP_STRUC.fitness, ...
        POP_STRUC.generation,ORG_STRUC.numGenerations,ORG_STRUC.stopCrit,...
        ORG_STRUC.stopFitness, ORG_STRUC.fitLimit, ORG_STRUC.optType);
    end
    if ToStop == 1
        fprintf(fp,'%30s \n', Message);
        Finish();
    else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp,'\n');
        WriteGenerationBackup();
        RATIO = POP_STRUC.good_and_fit / length(POP_STRUC.POPULATION);
        update_STUFF('INPUT.txt', RATIO, POP_STRUC.ranking);
        CreateCalcFolder();

        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); % same structure
        OFF_STRUC.POPULATION(1) = QuickStart(OFF_STRUC.POPULATION);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Operation = {'Heredity_M300', 'Random_M300', ...
             'Permutation_M300', 'SoftModeMutation_M300'};
        Num_Opera = [ORG_STRUC.howManyOffsprings, ORG_STRUC.howManyRand,...
            ORG_STRUC.howManyPermutations,  ORG_STRUC.howManyAtomMutations];
        count = 0;
        for i = 1 : length(Num_Opera)
            for j = 1 : Num_Opera(i)
                count = count + 1;
                eval([Operation{i} '(' num2str(count) ')']);
            end
        end

        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC.generation = POP_STRUC.generation + 1;

        fprintf(fp, [alignLine('Variation operators applied') '\n']);

        table = [strread(num2str(Num_Opera),'%s')'; Operation];
        fprintf(fp,' %10s structures produced by %s \n', table{:});

        addon_diff = KeepBestStructures();
        fprintf(fp,' %10d structures kept from previous gen\n',addon_diff);

        POP_STRUC = OFF_STRUC;
        num = pickUpSeeds();
        fprintf(fp, ' %10d structures added from Seeds/POSCARS\n',  num);
        fprintf(fp, [alignLine('Proceeding to the new generation ') '\n']);
        fprintf(fp, ' %10d parallel calculations running simutaneously\n',...
                     ORG_STRUC.numParallelCalcs);

        fprintf(fp, [alignLine( ...
                sprintf('Generation %d', POP_STRUC.generation) ) '\n'] );

        fprintf(fp,'  ID   Origin      Composition  Enthalpy(eV)');
        fprintf(fp,'  Volume(A^3)  KPOINTS  SYMMETRY\n');
        WriteGenerationStart();
        Start_POP();
        ORG_STRUC.currentGenDone = 0;  %avoid fitness calculation

        safesave ('Current_POP.mat', POP_STRUC)
        safesave ('Current_ORG.mat', ORG_STRUC)
        disp(' ');
        disp('New generation produced');

    end
end

fclose(fp);
