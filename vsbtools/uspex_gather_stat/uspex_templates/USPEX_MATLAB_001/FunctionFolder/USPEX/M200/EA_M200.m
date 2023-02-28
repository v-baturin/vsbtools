function EA_M200()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global USPEX_STRUC
global ANTISEEDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/OUTPUT.txt'];
fp = fopen(fpath, 'a+');
while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5

    if ORG_STRUC.currentGenDone == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LocalRelaxation();
        disp('status = Local optimisation finished')
        disp(' ')
        fprintf(fp, [alignLine('Local optimization finished') '\n']);
        fprintf(fp, '\n');
        fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d',...
                     POP_STRUC.generation) ) '\n'] );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 3:  Fitness                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SuccessRate();  %count how many structures were sucessfully done

        if ~((ORG_STRUC.pickedUP == 1) && ...
           (ORG_STRUC.pickUpGen == POP_STRUC.generation))
            fitness = CalcFitness_M200();
            Correlation(fitness);
            fprintf(fp, [alignLine(sprintf(...
            'Correlation coefficient = %.4f',...
            POP_STRUC.correlation_coefficient) ) '\n'] );

            % AntiSeedsCorrection(fitness);
            fitness = FitnessRankingCorrection(fitness);
            WriteGenerationOutput_M200(fitness);
            ORG_STRUC.currentGenDone=1;  %This is to avoid calculate fitness again
            safesave ([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC)
            safesave ([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC)
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ToStop = 0;
    if ~((ORG_STRUC.pickedUP == 1) && (ORG_STRUC.pickUpGen == POP_STRUC.generation))
        [ToStop, Message] = StopRun(POP_STRUC.fitness, POP_STRUC.generation, ...
        ORG_STRUC.numGenerations, ORG_STRUC.stopCrit, ORG_STRUC.stopFitness, ...
        ORG_STRUC.fitLimit, ORG_STRUC.optType);
    end
    if ToStop == 1
        fprintf(fp,'%30s \n', Message);
        Finish();
    else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp,'\n');

        WriteGenerationBackup();
        fraction_content = update_STUFF('INPUT.txt');
        fprintf(fp, [alignLine('VARIATION OPERATORS') '\n']);
        fprintf(fp, '    Keep %2d percent to produce next generation\n',...
                                             ORG_STRUC.bestFrac*100);
        fprintf(fp, fraction_content);
        fprintf(fp,'\n');

        CreateCalcFolder();

        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %make same structure
        OFF_STRUC.POPULATION(1) = QuickStart(OFF_STRUC.POPULATION);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Operation = {'Heredity_M200', 'Random_M200', 'Permutation_M200', ...
                      'SoftModeMutation_M200', 'LatticeMutation_M200'};
        Num_Opera = [ORG_STRUC.fracGene,ORG_STRUC.fracRand,ORG_STRUC.fracPerm,...
                     ORG_STRUC.fracAtomsMut, ORG_STRUC.fracLatMut];
        Num_Opera = round(ORG_STRUC.populationSize * Num_Opera);
 
        count = 0;
        for i = 1 : length(Num_Opera)
            for j = 1 : Num_Opera(i)
                count = count + 1;
                eval([Operation{i} '(' num2str(count) ')']);
            end
        end

        if exist('spinOperation')
            spinOperation(count);
        end

        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC.generation = POP_STRUC.generation + 1;

        fprintf(fp, [alignLine('Variation operators applied') '\n']);

        table = [strread(num2str(Num_Opera), '%s')'; Operation];
        fprintf(fp,' %10s structures produced by %s \n', table{:});

        addon_diff = KeepBestStructures();
        fprintf(fp,' %10d structures kept from previous generation\n', addon_diff);

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
        ORG_STRUC.currentGenDone = 0;  %important to avoid fitness calculation

        safesave ('Current_POP.mat', POP_STRUC)
        safesave ('Current_ORG.mat', ORG_STRUC)
        disp(' ');
        disp('New generation produced');

    end
end

fclose(fp);
