function MINHOP()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC
global USPEX_STRUC

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%     Local Relaxation     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MINHOP_LocalRelaxation();
    [ranking, bad_rank, fitness] = MINHOP_CalcFitness(ORG_STRUC.conv_till);
    %  MINHOP_WriteGenerationOutput(fitness, 1);
    
    if ~((ORG_STRUC.pickedUP == 1) && (ORG_STRUC.pickUpGen == POP_STRUC.generation))
        if (ORG_STRUC.FullRelax > 0) && (POP_STRUC.finalOptimization == 0)
            fprintf(fp, [alignLine('After Relaxation') '\n']);

            if ORG_STRUC.varcomp == 1
                fprintf(fp,'  ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY\n');
            else
                fprintf(fp,'  ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY\n');
            end
            
            POP_STRUC.finalOptimization = 1;
            for i = 1 : length(POP_STRUC.POPULATION) %- bad_rank   mahdi
                POP_STRUC.POPULATION(ranking(i)).Done = 0;
                POP_STRUC.POPULATION(ranking(i)).ToDo = 1;
            end
            MINHOP_LocalRelaxation();
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%     Local Relaxation     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fp, [alignLine('Local optimization finished') '\n']);
    fprintf(fp, ['\n']);
    fprintf(fp, ['\n']);
    total_population = length(USPEX_STRUC.POPULATION);
    if total_population > 2
        the_same_structure = 2;
    elseif total_population < 2
        the_same_structure = 0;
    end
    
    %%%%%%To check the structure with all found structures with finger print criteria and energy difference check!
    %                 fingerprint001 = USPEX_STRUC.POPULATION(total_population).FINGERPRINT;
    fingerprint001 = POP_STRUC.POPULATION(1).FINGERPRINT;
    E1 = POP_STRUC.POPULATION(1).Enthalpies(end)/sum(POP_STRUC.POPULATION(1).numIons);
    cosdis = zeros(total_population,1);
    for a = 1 : total_population -1
        E2 = USPEX_STRUC.POPULATION(a).Enthalpies(end)/sum(USPEX_STRUC.POPULATION(a).numIons);
        diff = abs(E1-E2);
        
        
        cosdis(a) = cosineDistance(fingerprint001, USPEX_STRUC.POPULATION(a).FINGERPRINT, ORG_STRUC.weight);
        if  (cosdis(a) < ORG_STRUC.toleranceFing) &&  (diff < 0.015)
            the_same_structure = 1;
            disp('The same after relaxation!')
            
            break;
        else
            the_same_structure = 0;
            
        end
    end
    
    if  the_same_structure == 0
        MINHOP_WriteGenerationOutput(fitness, 0);
        MINHOP_GenerationBackup();
    end
    ToStop = 0;
    if POP_STRUC.generation >= ORG_STRUC.numGenerations
        ToStop = 1;
    end
    
    if ToStop==1
        MINHOP_Finish();
    else
        
        ind = ranking(1);
        COORDINATES    = POP_STRUC.POPULATION(ind).coords0;
        LATTICE        = POP_STRUC.POPULATION(ind).lat0;
        PressureTensor = POP_STRUC.POPULATION(ind).PressureTensor0;
        
        % Reduce supercell, if needed. We skip it in the first generation.
        if ORG_STRUC.varcomp
            [reduced, newCoord, newLat, newNumIons, newSupercell] = reduceSupercell(LATTICE, COORDINATES, ...
                POP_STRUC.POPULATION(ind).numIons, POP_STRUC.POPULATION(ind).superCell, 0.3);
            if reduced
                fprintf(fp, 'Translational symmetry detected. The supercell size was reduced.\n');
                POP_STRUC.POPULATION(ind).LATTICE = newLat;
                POP_STRUC.POPULATION(ind).COORDINATES = newCoord;
                POP_STRUC.POPULATION(ind).numIons = newNumIons;
                POP_STRUC.POPULATION(ind).superCell = newSupercell;
            end
        end
        CreateCalcFolder();
        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %just make it has same structure
        OFF_STRUC.POPULATION(1) = POP_STRUC.POPULATION(ind);
        OFF_STRUC.POPULATION(1).howCome = 'keptBest';
        OFF_STRUC.POPULATION(1).Parents.parent = num2str(ind);
        OFF_STRUC.POPULATION(1).Parents.enthalpy = 0;
        OFF_STRUC.POPULATION(1).Step = [];
        
        toSkip = 0;
        for i = 2 : ORG_STRUC.populationSize
            if toSkip
                toSkip = 0;
                continue
            end
            
            if ORG_STRUC.varcomp
                disp('NOT DONE YET')
                %atomSoftModeMutation_varcomp(i, ind, 0)
            else
                
                MINHOP_SoftModeMutation(i, ind, the_same_structure)
                
            end
            if length(OFF_STRUC.POPULATION) > ORG_STRUC.populationSize
                OFF_STRUC.POPULATION(i+1) = OFF_STRUC.POPULATION(end);
                OFF_STRUC.POPULATION(end) = [];
                toSkip = 1;
            end
        end
        
        if ORG_STRUC.varcomp
            %%%%%%%%%%%%%%  if the 'basic' cell is chosen as parent alternative:  %%%%%%%%%
            %%%%%%%%%%%%%%         we want to relax this basic cell first         %%%%%%%%%
            if (mod(POP_STRUC.generation, ORG_STRUC.useBasicCell) == ORG_STRUC.useBasicCell-1) && (sum(POP_STRUC.POPULATION(ind).superCell) > 3)
                OFF_STRUC.POPULATION(end+1) = POP_STRUC.POPULATION(ind);
                OFF_STRUC.POPULATION(end).COORDINATES = ORG_STRUC.bestBasicStructure;
                OFF_STRUC.POPULATION(end).Step = [];
                OFF_STRUC.POPULATION(end).numIons = ORG_STRUC.numIons;
                OFF_STRUC.POPULATION(end).Softmode_num = 0;
                POP_STRUC.basicStructureNumber = length(OFF_STRUC.POPULATION); % remember the basic structure that we relax
            elseif (mod(POP_STRUC.generation, ORG_STRUC.useBasicCell) == 0) && (sum(POP_STRUC.POPULATION(ind).superCell) > 3) && (ind < POP_STRUC.basicStructureNumber)
                disp(' ');
                disp('Softmutating the basic cell, to increase diversity');
                disp(' ');
                populSize = length(OFF_STRUC.POPULATION);
                OFF_STRUC.POPULATION(end+1) = POP_STRUC.POPULATION(ind);
                if POP_STRUC.DoneOrder(POP_STRUC.basicStructureNumber) > 0
                    OFF_STRUC.POPULATION(end).COORDINATES = POP_STRUC.POPULATION(POP_STRUC.basicStructureNumber).COORDINATES;  % relaxed basic structure
                else
                    OFF_STRUC.POPULATION(end).COORDINATES = ORG_STRUC.bestBasicStructure;  % old algorithm
                end
                OFF_STRUC.POPULATION(end).numIons = ORG_STRUC.numIons;
                OFF_STRUC.POPULATION(populSize+1).superCell = [1 1 1];
                POP_STRUC.POPULATION(end+1) = OFF_STRUC.POPULATION(populSize + 1);
                
                toSkip = 0;
                for i = 2 : ORG_STRUC.populationSize
                    if toSkip
                        toSkip = 0;
                        continue
                    end
                    disp('NOT DONE YET')
                    
                    %          atomSoftModeMutation_varcomp(i+populSize, 1+populSize, 1)
                    % fix the number of structures, if both softmutation directions were used
                    if length(OFF_STRUC.POPULATION) > populSize+ORG_STRUC.populationSize
                        OFF_STRUC.POPULATION(i+ORG_STRUC.populationSize+1) = OFF_STRUC.POPULATION(end);
                        OFF_STRUC.POPULATION(end) = [];
                        toSkip = 1;
                    end
                end
            end
        end
        
        i = 1;
        while i <= length(OFF_STRUC.POPULATION)
            if isempty(OFF_STRUC.POPULATION(i).superCell)
                OFF_STRUC.POPULATION(i) = [];
            else
                i = i + 1;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         OFF_STRUC.generation = POP_STRUC.generation + 1;
         POP_STRUC = OFF_STRUC;
         if the_same_structure ~= 1
          WriteGenerationStart();
         end
         MINHOP_Start_POP();

        fprintf(fp, [alignLine( sprintf('Generation %d', POP_STRUC.generation) ) '\n'] );
        fprintf(fp, [alignLine('Before Relaxation') '\n']);

        if ORG_STRUC.varcomp == 1
            fprintf(fp,'  ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY\n');
        else
            fprintf(fp,'  ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY\n');
        end
        
        save ('Current_POP.mat', 'POP_STRUC')
        save ('Current_ORG.mat', 'ORG_STRUC')
        disp('New generation produced');
    end
end
fclose(fp);
