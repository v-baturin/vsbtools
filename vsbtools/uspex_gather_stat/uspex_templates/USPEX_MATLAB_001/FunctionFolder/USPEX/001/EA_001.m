function EA_001()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global USPEX_STRUC
global ANTISEEDS
global CLUSTERS
global AR_VARIANTS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
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
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', ...
                                    POP_STRUC.generation) ) '\n'] );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 3:  Fitness                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SuccessRate();  %count how many structures were sucessfully done
    
        if ~((ORG_STRUC.pickedUP == 1) && ...
            (ORG_STRUC.pickUpGen == POP_STRUC.generation))
            fitness = CalcFitness_001();
            Correlation(fitness);
            fprintf(fp, [alignLine(sprintf(...
            'Correlation coefficient = %.4f',...
            POP_STRUC.correlation_coefficient) ) '\n'] );
            
            fitness = AntiSeedsCorrection_001(fitness);
            fitness = FitnessRankingCorrection(fitness);
            WriteClusters(fitness);
            WriteGenerationOutput_001(fitness);
        end
        ORG_STRUC.currentGenDone=1;  %This is to avoid we calculate fitness again
        safesave ([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC)
        safesave ([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ToStop = 0;
    if ~((ORG_STRUC.pickedUP == 1) && (ORG_STRUC.pickUpGen == POP_STRUC.generation))
        [ToStop, Message] = StopRun(POP_STRUC.fitness, POP_STRUC.generation, ...
        ORG_STRUC.numGenerations, ORG_STRUC.stopCrit, ORG_STRUC.stopFitness, ...
        ORG_STRUC.fitLimit, ORG_STRUC.opt_sign * ORG_STRUC.optType);
    end
    if ToStop == 1
        fprintf(fp,'%30s \n', Message);
        Finish();
    else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
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
        InitializeVariantsAR();
        plot_all_isomers();
        
        Operation = {'Heredity_001', 'H1eredity_001', 'Random_001', ...
                  'Permutation_001', 'Transmutation_001', 'SoftModeMutation_001', ...
                  'AddAtom_001', 'RemoveAtom_001'};
        Num_Opera = [ORG_STRUC.fracGene, ORG_STRUC.fracG1ene, ORG_STRUC.fracRand, ...
                 ORG_STRUC.fracPerm, ORG_STRUC.fracTrans, ORG_STRUC.fracAtomsMut, ...
                 ORG_STRUC.fracAddAtom, ORG_STRUC.fracRemAtom];
        Num_Opera = round(ORG_STRUC.populationSize * Num_Opera);
        count = 0;
        for i = 1 : length(Num_Opera)
            for j = 1:Num_Opera(i)
                count = count + 1;
                eval([Operation{i} '(' num2str(count) ')']);
            end
        end
        safesave ('CLUSTERS.mat', CLUSTERS)
        
        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC.generation = POP_STRUC.generation + 1;
        
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Variation operators applied') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);

        table = [strread(num2str(Num_Opera),'%s')'; Operation];
        fprintf(fp,' %10s structures produced by %s \n', table{:});

        addon_diff = KeepBestStructures_001();
        fprintf(fp,' %10d structures kept from previous generation\n',addon_diff);
        
        POP_STRUC = OFF_STRUC;
        num = pickUpSeeds();
        fprintf(fp,' %10d structures added from Seeds/POSCARS\n',  num);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Proceeding to the new generation ') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp,' %10d parallel calculations running simutaneously\n',...
                                             ORG_STRUC.numParallelCalcs);
        fprintf(fp, [alignLine('-', 0) '\n']);
        
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
        
        %%% delete n-5 generationa folder in results
        if POP_STRUC.generation > 5
            unix(['rm -rf ' ORG_STRUC.resFolder '/generation' num2str(POP_STRUC.generation-5)]);
        end
    end
end

fclose(fp);       
