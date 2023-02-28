function PSO()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global ANTISEEDS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     Evolutionary Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 2:        Local Relaxation    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    while 1 % this eternal cycle is needed for non parallel calculations, parallelized one will /break out of it

        ReadJobs_PSO();
        SubmitJobs_PSO();

        if ORG_STRUC.platform > 0 %
            if sum([POP_STRUC.POPULATION(:).Done]) ~= length(POP_STRUC.POPULATION)
                delete('still_running');
                fclose ('all');
                quit
            else
                break;
            end
        else
            if sum([POP_STRUC.POPULATION(:).Done]) == length(POP_STRUC.POPULATION)
                break; % break out of eternal cycle when non parallel calculations finished for given generation
            end
        end
    end

    disp('status = Local optimization finished')
    disp(' ')
    fprintf(fp, [alignLine('-', 0) '\n']);
    fprintf(fp, [alignLine('Local optimization finished') '\n']);
    fprintf(fp, '\n');
    fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', POP_STRUC.generation) ) '\n'] );

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 3:  Fitness                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    good_and_fit = 0;
    for fit_loop = 1 : length(POP_STRUC.POPULATION)
        if POP_STRUC.POPULATION(fit_loop).Enthalpies(end) < 90000
            good_and_fit = good_and_fit + 1;
        end
    end

    if good_and_fit < floor(length(POP_STRUC.POPULATION) / 3)
        fprintf(fp, 'Too many structures have errors or failed the constraints after optimization.\n');
        fprintf(fp, 'Please check the input files. The calculation has to stop.\n');
        fprintf(fp, 'Possible reasons: badly tuned optimization parameters or unreasonable contraints.\n');
        quit;
    elseif good_and_fit / length(POP_STRUC.POPULATION) < ORG_STRUC.bestFrac
        fprintf(fp, 'Some structures have errors or failed the constraints after optimisation,\n');
        fprintf(fp, 'selection of these structures as parents is possible.\n');
        fprintf(fp, 'bestFrac parameter will be lowered for this generation to discard such structures.\n');
        fprintf(fp, 'Please check the input files. Results may not be reliable. \n');
        fprintf(fp, 'Possible reasons: high bestFrac, bad optimization parameters or unreasonable contraints.\n');
    end

    if ~((ORG_STRUC.pickedUP == 1) & (ORG_STRUC.pickUpGen == POP_STRUC.generation))
        fitness = CalcFitness_PSO();
        cor_coefficient = Correlation(fitness);
        fprintf(fp, [alignLine( sprintf('Correlation coefficient = %.4f', cor_coefficient) ) '\n'] );
        
        AntiSeedsCorrection(fitness);
        fitness = FitnessRankingCorrection(fitness);
        WriteGenerationOutput_PSO(fitness);

        for ind = 1 : ORG_STRUC.populationSize
            if (POP_STRUC.generation == 1) || (POP_STRUC.PSO(ind).fitness > fitness(ind))
                POP_STRUC.PSO(ind).lattice = POP_STRUC.POPULATION(ind).LATTICE;
                POP_STRUC.PSO(ind).coordinates = POP_STRUC.POPULATION(ind).COORDINATES;
                POP_STRUC.PSO(ind).order = POP_STRUC.POPULATION(ind).order;
                POP_STRUC.PSO(ind).numIons = POP_STRUC.POPULATION(ind).numIons;
                POP_STRUC.PSO(ind).fitness = fitness(ind);
                POP_STRUC.PSO(ind).fingerprint = POP_STRUC.POPULATION(ind).FINGERPRINT;
                POP_STRUC.PSO(ind).enthalpy = POP_STRUC.POPULATION(ind).Enthalpies(end);
                POP_STRUC.PSO(ind).Number = POP_STRUC.POPULATION(ind).Number;
            end
            if POP_STRUC.PSO(POP_STRUC.bestPSOstruc).fitness > POP_STRUC.PSO(ind).fitness
                POP_STRUC.bestPSOstruc = ind;
            end
        end
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 4:  If stop this calculation                 %%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ToStop=0;
    if POP_STRUC.generation >= ORG_STRUC.numGenerations
        ToStop = 1;
    elseif POP_STRUC.generation >= ORG_STRUC.stopCrit
        same_in = SameBest(fitness);
        if same_in == ORG_STRUC.stopCrit - 1
            fprintf(fp,'Halting criteria achieved: %d\n', POP_STRUC.bodyCount);
            ToStop = 1;
        end
    elseif ~isempty(ORG_STRUC.stopFitness)
        if fitness(ranking(1)) <= ORG_STRUC.stopFitness
            fprintf(fp,'Known ground state achieved: %d\n', POP_STRUC.bodyCount);
            ToStop = 1;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Step 5: Selection and Variation               %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ToStop==1
        Finish_PSO();
    else
        fprintf(fp, [alignLine('Proceeding to Selection Process') '\n']);
        fprintf(fp,'\n');

        WriteGenerationBackup();

        update_STUFF('INPUT.txt');
        CreateCalcFolder();

        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %just make it has same structure
        OFF_STRUC.POPULATION(1) = QuickStart(OFF_STRUC.POPULATION);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ind = 1 : ORG_STRUC.populationSize
            dist1 = cosineDistance(POP_STRUC.POPULATION(ind).FINGERPRINT, POP_STRUC.PSO(ind).fingerprint, ORG_STRUC.weight);   %Dp in 2013 CPC paper
            dist2 = cosineDistance(POP_STRUC.POPULATION(ind).FINGERPRINT, POP_STRUC.PSO(POP_STRUC.bestPSOstruc).fingerprint, ORG_STRUC.weight); %Dg in 2013 CPC paper
            P_p = rand * ORG_STRUC.PSO_BestEver * dist1;
            P_g = rand * ORG_STRUC.PSO_BestStruc * dist2;
            P_m = rand * ORG_STRUC.PSO_softMut;
            P_Normal = (P_p + P_g + P_m);
            P_r = ORG_STRUC.fracRand;
            if rand < P_r | POP_STRUC.POPULATION(ind).Error > ORG_STRUC.maxErrors %first random number is used to check for random generation
                Random_PSO(ind);
            else
                if rand < P_m/P_Normal
                    SoftModeMutation_PSO(ind);
                else
                    Heredity_PSO(ind); % there we randomly choose between 2 heredity types
                end
            end
        end

        OFF_STRUC.SOFTMODEParents = POP_STRUC.SOFTMODEParents;
        disp(' ');
        disp('Variation operators applied, applying elitist scheme')

        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC.generation = POP_STRUC.generation + 1;

        fprintf(fp, [alignLine('Variation operators applied') '\n']);

        addon_diff = KeepBestStructures_PSO();
        fprintf(fp,'            %4d structures kept as best from the previous generation\n', addon_diff);

        POP_STRUC = OFF_STRUC;
        num = pick_Seeds();
        fprintf(fp,'            %4d Seeds structures are added from Seeds/POSCARS\n', num);
        fprintf(fp, [alignLine('Proceeding to the new generation relaxation') '\n']);
        fprintf(fp,'            %4d parallel calculations are performed simutaneously\n', ORG_STRUC.numParallelCalcs);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf(fp, [alignLine( sprintf('Generation %d', POP_STRUC.generation) ) '\n'] );
        fprintf(fp,'  ID   Origin      Composition  Enthalpy(eV)  Volume(A^3)  KPOINTS  SYMMETRY\n');
        WriteGenerationStart();
        Start_POP_PSO();

        save ('Current_POP.mat', 'POP_STRUC')
        save ('Current_ORG.mat', 'ORG_STRUC')
        disp(' ');
        disp('New generation produced');
    end
end
fclose(fp);
