function EA_301()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global POOL_STRUC
global ANTISEEDS
global USPEX_STRUC


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');
while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5
    if ORG_STRUC.currentGenDone == 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LocalRelaxation();
        disp('status = Local optimisation finished')
        disp(' ')
        fprintf(fp, [alignLine('Local optimization finished') '\n']);
        fprintf(fp, '\n');
        fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', ...
                    POP_STRUC.generation) ) '\n'] );

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 3:  Fitness                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SuccessRate();  %count how many structures were sucessfully done

        if ~((ORG_STRUC.pickedUP == 1) && ...
           (ORG_STRUC.pickUpGen == POP_STRUC.generation))
           fitness = CalcFitness_301();
           Correlation(fitness);
           fprintf(fp, [alignLine(sprintf(...
           'Correlation coefficient = %.4f',...
           POP_STRUC.correlation_coefficient) ) '\n'] );

           fitness = AntiSeedsCorrection(fitness);
           if ORG_STRUC.paretoRanking ~= 0
              [ConVeX, convex_hull0] = update_convex_hull(ORG_STRUC.numIons);
	      for a = 1 : length(POP_STRUC.POPULATION) 
	          ConVeX(a) = ConVeX(a) * sum(POP_STRUC.POPULATION(a).numBlocks) / sum(POP_STRUC.POPULATION(a).numIons); % eV/atom
	      end
              fitness = ParetoFitnessRankingCorrection(fitness, ConVeX);
           else
              fitness = FitnessRankingCorrection(fitness);
           end
           POOL_STRUC.Composition_ranking = ...
                          zeros(length(POOL_STRUC.Composition_ranking),1);
           POOL_STRUC.POPULATION  = [];
           POOL_STRUC.POPULATION  = struct('COORDINATES', {}, ...
            'LATTICE', {},  'numIons', {}, 'numBlocks',{},  ...
            'order', {}, 'enthalpy', {}, 'Number', {}, 'FINGERPRINT', {});

           compositionStatistic(ORG_STRUC.firstGeneSplitAll);
           WriteGenerationOutput_301(fitness);
        end
        WriteGenerationBackup();
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Step 4:  Update the USPEX.mat and POOL.mat with AntiSeeds     %%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % reupdateStructures_301();

        %Step 5: save and quit

        ORG_STRUC.currentGenDone=1;
        safesave ([ORG_STRUC.homePath '/Current_POP.mat'],  POP_STRUC)
        safesave ([ORG_STRUC.homePath '/Current_ORG.mat'],  ORG_STRUC)

        if POP_STRUC.generation >= ORG_STRUC.numGenerations
            Finish();
            fclose('all');
            return;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 5:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ToStop = 0;
    if ~((ORG_STRUC.pickedUP == 1) && (ORG_STRUC.pickUpGen == POP_STRUC.generation))
      [ToStop, Message] = StopRun(POP_STRUC.fitness, POP_STRUC.generation,...
      ORG_STRUC.numGenerations, ORG_STRUC.stopCrit, ORG_STRUC.stopFitness,...
      ORG_STRUC.fitLimit, ORG_STRUC.opt_sign * ORG_STRUC.optType);
    end
    if ToStop == 1
        fprintf(fp,'%30s \n', Message);
        Finish();
    else

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 6: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Start the new generation %%
        waitForCOPEX();
        
        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp, '\n');
        if ~isempty(POP_STRUC.convex_hull)
            total_comp = sum(POOL_STRUC.Composition_ranking(:) > 0);
            total_pop = length(POOL_STRUC.POPULATION);
        
            fprintf(fp,'N_compositions in selection pool: %4d\n', total_comp);
            fprintf(fp,'N_structures in selection pool: %4d\n\n', total_pop);
            fprintf(fp,'Compositions      N_Structures     BestEnthalpies(eV/atom)  fitness    Surviving Gen\n');
            for i=1:length(POOL_STRUC.Composition_ranking)
                comp = num2str(POOL_STRUC.Composition_ratio(i, :), '%6.3f');
                rank = POOL_STRUC.Composition_ranking(i);
                best = POOL_STRUC.Composition_Bestenthalpy(i);
                surv = POOL_STRUC.Composition_surviving(i);
                fit  = POOL_STRUC.Composition_fitness(i);
                if rank > 0
                    fprintf(fp, '%-18s    %4d         %10.4f     %16.3f       %4d\n', comp, rank, best, fit, surv);
                end
            end
        end
        
        fprintf(fp,'\n');
        fraction_content = update_STUFF('INPUT.txt');
        fprintf(fp, [alignLine('VARIATION OPERATORS') '\n']);
        fprintf(fp, '    Keep %2d percent to produce next generation\n',...
                                             ORG_STRUC.bestFrac*100);
        fprintf(fp, fraction_content);
        fprintf(fp, '\n');

        if size(ORG_STRUC.atomType, 2) < 4  %for less than 4 component
            coordinating_number_analisis(); %determine coordinating numbers
        end                                 %in topological random

        poolsize = length(POOL_STRUC.POPULATION);

        ORG_STRUC.tournament = zeros(poolsize, 1);
        ORG_STRUC.tournament(poolsize) = 1;
        for loop = 2 : poolsize
            ORG_STRUC.tournament(end-loop+1) = ...
            ORG_STRUC.tournament(end-loop+2) + loop ^ 2;
        end

        if POP_STRUC.generation == 1
            if poolsize > ORG_STRUC.initialPopSize
                ORG_STRUC.tournament = zeros(ORG_STRUC.initialPopSize, 1);
                ORG_STRUC.tournament(ORG_STRUC.initialPopSize) = 1;
                for loop = 2 : ORG_STRUC.initialPopSize
                    ORG_STRUC.tournament(end-loop+1) = ...
                    ORG_STRUC.tournament(end-loop+2) + loop ^ 2;
                end
            end
        end
       
        %%% This part was added for pareto, to give all the structures
        %%% of same Pareto front, the same posibility to get chosen!
        if ORG_STRUC.paretoRanking ~= 0
            TotalStructures = sum(POOL_STRUC.paretoFront);  %% equal to length(POOL_STRUC.POPULATION)
            Tournament      = zeros(TotalStructures , 1);
            for a = 1 : length(ORG_STRUC.tournament)
                Tournament(a) = ORG_STRUC.tournament(a);
            end
            EndPoint = 0;
            for b = 1 : length(POOL_STRUC.paretoFront)
                StartPiont = 2 + EndPoint;
                EndPoint   = sum(POOL_STRUC.paretoFront(1:b));
                Tournament(StartPiont : EndPoint) = 0;
            end
            for a = length(Tournament) : -1 : 1
                if Tournament(a) == 0
                    Tournament(a) = [];
                end
            end
            ORG_STRUC.tournament = Tournament;
        end
        CreateCalcFolder();
        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %make same struct
        OFF_STRUC.POPULATION(1) = QuickStart(OFF_STRUC.POPULATION);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Operation = {'Heredity_301', 'Random_301', 'RandTop_301', ...
                     'Permutation_301', 'Transmutation_301', ...
                     'SoftModeMutation_301', 'LatticeMutation_301'};
        Num_Opera = [ORG_STRUC.fracGene, ORG_STRUC.fracRand,...
               ORG_STRUC.fracRandTop,  ORG_STRUC.fracPerm, ...
               ORG_STRUC.fracTrans, ORG_STRUC.fracAtomsMut, ...
                                      ORG_STRUC.fracLatMut];
        Num_Opera = round(ORG_STRUC.populationSize * Num_Opera);
        count = 1;
        for i = 1 : length(Num_Opera)
            for j = 1 : Num_Opera(i)               
                eval([Operation{i} '(' num2str(count) ')']);
                count = length(OFF_STRUC.POPULATION) + 1;
            end
        end
        
        count = length(OFF_STRUC.POPULATION) + 1;
        if exist('spinOperation')
            spinOperation(count);
        end
        
        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme')
        
        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC.generation = POP_STRUC.generation + 1;
        
        fprintf(fp, [alignLine('Variation operators applied') '\n']);
        table = [strread(num2str(Num_Opera),'%s')'; Operation];
        fprintf(fp,' %10s structures produced by %s \n', table{:});
        
        fprintf(fp,' %10d structures produced by spinmutation\n',...
                    ORG_STRUC.howManySpinmutations);
        
        addon_diff = KeepBestStructures();
        fprintf(fp,' %10d structures kept from previous gen\n',  addon_diff);
        [addon_copex,OFF_STRUC] = importCOPEXStructures(...
                                  OFF_STRUC, POP_STRUC.generation);
        fprintf(fp,' %10d structures from other USPEX Calcs\n', addon_copex);
        
        POP_STRUC = OFF_STRUC;
%       num = pick_Seeds();
        num = pickUpSeeds();
        fprintf(fp,' %10d structures added from Seeds/POSCARS\n', num);
        fprintf(fp, [alignLine('Proceeding to a new generation') '\n']);
        fprintf(fp,' %10d parallel calculations running simultaneously\n',...
                                                ORG_STRUC.numParallelCalcs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf(fp,'  ID   Origin      Composition  Enthalpy(eV)');
        fprintf(fp,'  Volume(A^3)  KPOINTS  SYMMETRY\n');
        WriteGenerationStart();
        Start_POP();
        ORG_STRUC.currentGenDone = 0;  %important to avoid fitness calculation
        
        if ( ORG_STRUC.pluginType > 0 ) && ...
           ( mod( POP_STRUC.generation, ORG_STRUC.pluginType)==0 )
            ORG_STRUC.startNextGen = 0;
        else
            ORG_STRUC.startNextGen = 1;
        end
        
        safesave ('Current_POP.mat', POP_STRUC)
        safesave ('Current_ORG.mat', ORG_STRUC)
        
        disp(' ');
        disp('New generation produced');
        %Start the new generation
    end   %while generation cycle
end

fclose(fp);

