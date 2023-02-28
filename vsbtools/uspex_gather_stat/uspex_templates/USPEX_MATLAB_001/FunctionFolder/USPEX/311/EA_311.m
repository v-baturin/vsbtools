function EA_311()

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
global POOL_STRUC
global OFF_STRUC
global ANTISEEDS

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
        fprintf(fp, [alignLine('-', 0) '\n']);
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
            fitness = CalcFitness_311();
            Correlation(fitness);
            fprintf(fp, [alignLine(sprintf(...
            'Correlation coefficient = %.4f',...
            POP_STRUC.correlation_coefficient) ) '\n'] );
            
            AntiSeedsCorrection(fitness);
            fitness = FitnessRankingCorrection(fitness);
            
            POOL_STRUC.Composition_ranking = ...
                               zeros(length(POOL_STRUC.Composition_ranking),1);
            POOL_STRUC.POPULATION  = [];
            POOL_STRUC.POPULATION  = struct('COORDINATES', {}, 'LATTICE', {},...
                                'numIons', {}, 'numBlocks',{}, 'numMols', {},...
                                   'order', {}, 'enthalpy', {}, 'Number', {});
            
            compositionStatistic(ORG_STRUC.firstGeneSplitAll);
            WriteGenerationOutput_311(fitness);
        end
        ORG_STRUC.currentGenDone=1;  %This is to avoid we    
        safesave ([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC)
        safesave ([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC)
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if POP_STRUC.generation >= ORG_STRUC.numGenerations
        Finish();
    else
        
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp, [alignLine('-', 0) '\n']);
        fprintf(fp,'\n');
        if ~isempty(POP_STRUC.convex_hull)
            total_comp = sum(POOL_STRUC.Composition_ranking(:)>0);
            fprintf(fp,'The number of compositions in the current pool: %4d\n',...
                                                                   total_comp);
            total_pop = length(POOL_STRUC.POPULATION);
            fprintf(fp,'The number of structures in the current pool: %4d\n\n',...
                                                                    total_pop);
            fprintf(fp,'Compositions      N_Structures     BestEnthalpies');
            fprintf(fp,'(eV/atom)  fitness    Surviving Gen\n');
            for i=1:length(POOL_STRUC.Composition_ranking)
                comp = num2str(POOL_STRUC.Composition_ratio(i,:), '%6.3f');
                rank = POOL_STRUC.Composition_ranking(i);
                best = POOL_STRUC.Composition_Bestenthalpy(i);
                surv = POOL_STRUC.Composition_surviving(i);
                fit = POOL_STRUC.Composition_fitness(i);
                if rank > 0
                    fprintf(fp,'%-18s %8d %18.4f %20.3f %10d\n',...
                                comp, rank, best, fit, surv);
                end
            end
        end
        
        POOL2MOL();
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        WriteGenerationBackup();
        fraction_content = update_STUFF('INPUT.txt');
        fprintf(fp, [alignLine('VARIATION OPERATORS') '\n']);
        fprintf(fp, '    Keep %2d percent to produce next generation\n',...
                                             ORG_STRUC.bestFrac*100);
        fprintf(fp, fraction_content);
        fprintf(fp,'\n');

        poolsize = length(POOL_STRUC.POPULATION);
        ORG_STRUC.tournament = zeros(poolsize,1);
        ORG_STRUC.tournament(poolsize) = 1;
        for loop = 2:poolsize
            ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2) + loop^2;
        end
        
        if POP_STRUC.generation == 1
            if poolsize > ORG_STRUC.initialPopSize
                ORG_STRUC.tournament = zeros(ORG_STRUC.initialPopSize,1);
                ORG_STRUC.tournament(ORG_STRUC.initialPopSize) = 1;
                for loop = 2:ORG_STRUC.initialPopSize
                    ORG_STRUC.tournament(end-loop+1) = ORG_STRUC.tournament(end-loop+2)+loop^2;
                end
            end
        end
        
        CreateCalcFolder();
        
        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %make same structure
        OFF_STRUC.POPULATION(1) = QuickStart(OFF_STRUC.POPULATION);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Operation = {'Heredity_311', 'Random_311', 'Permutation_311', ...
            'SoftModeMutation_311', 'LatMutation_311', 'Rotation_311'};
        Num_Opera = [ORG_STRUC.fracGene,ORG_STRUC.fracRand,ORG_STRUC.fracPerm, ...
            ORG_STRUC.fracAtomsMut, ORG_STRUC.fracLatMut, ORG_STRUC.fracRotMut ];
        Num_Opera = round(ORG_STRUC.populationSize * Num_Opera);
        count = 0;
        for i = 1 : length(Num_Opera)
            for j = 1:Num_Opera(i)
                count = count + 1;
                eval([Operation{i} '(' num2str(count) ')']);
            end
        end
        
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

        addon_diff = KeepBestStructures();
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
        
    end
end

fclose(fp);           
