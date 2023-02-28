function META()

global ORG_STRUC
global POP_STRUC
global OFF_STRUC

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%     THE ACTUAL ALGORITHM    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while POP_STRUC.generation < ORG_STRUC.numGenerations + 0.5

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% 1,     Local Relaxation     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    META_LocalRelaxation();
    [ranking, bad_rank, fitness] = META_CalcFitness(ORG_STRUC.conv_till);
    if POP_STRUC.finalOptimization == 0
        META_WriteGenerationOutput(fitness, 1);
    end
    %%%Switch to Full Relax if needs,
    if ~((ORG_STRUC.pickedUP == 1) && (ORG_STRUC.pickUpGen == POP_STRUC.generation))
        if (ORG_STRUC.FullRelax > 0) && (POP_STRUC.finalOptimization == 0)
            fprintf(fp, [alignLine('Full Relax') '\n']);

            if ORG_STRUC.varcomp == 1
                fprintf(fp,'  ID    SuperCell   Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY\n');
            else
                fprintf(fp,'  ID   Composition  Enthalpy(eV/atom)  Volume(A^3/atom)  KPOINTS  SYMMETRY\n');
            end

            POP_STRUC.finalOptimization = 1;
            if ORG_STRUC.FullRelax == 1  % 'final' optimization only for best structure
                POP_STRUC.POPULATION(ranking(1)).Done = 0;
                POP_STRUC.POPULATION(ranking(1)).ToDo = 1;
            else
                for i = 1 : length(POP_STRUC.POPULATION) - bad_rank
                    POP_STRUC.POPULATION(ranking(i)).Done = 0;
                    POP_STRUC.POPULATION(ranking(i)).ToDo = 1;
                end
            end
            META_LocalRelaxation();
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% 2,     Local Relaxation     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp('status = Local optimisation finished')
    fprintf(fp, [alignLine('Local optimization finished') '\n']);
    fprintf(fp, '\n');
    fprintf(fp, [alignLine( sprintf('SUMMARY of Generation %d', POP_STRUC.generation) ) '\n'] );

    META_WriteGenerationOutput(fitness, 0);
    META_GenerationBackup();
    createORG_META('INPUT.txt');

    ToStop = 0;
    if POP_STRUC.generation >= ORG_STRUC.numGenerations
        ToStop = 1;
    end

    if ToStop == 1
        META_Finish();
    else
        ind = ranking(1);
        COORDINATES    = POP_STRUC.POPULATION(ind).coords0;
        LATTICE        = POP_STRUC.POPULATION(ind).lat0;
        PressureTensor = POP_STRUC.POPULATION(ind).PressureTensor0;

        % Reduce supercell, if needed. We skip it in the first generation.
        if ORG_STRUC.varcomp
            [reduced, newCoord, newLat, newNumIons, newSupercell] = reduceSupercell(LATTICE, COORDINATES, ...
                POP_STRUC.POPULATION(ind).numIons, POP_STRUC.POPULATION(ind).superCell, 0.2);
            if reduced
                fprintf(fp, 'Translational symmetry detected. The supercell size was reduced.\n');
                POP_STRUC.POPULATION(ind).LATTICE = newLat;
                POP_STRUC.POPULATION(ind).COORDINATES = newCoord;
                POP_STRUC.POPULATION(ind).numIons = newNumIons;
                POP_STRUC.POPULATION(ind).superCell = newSupercell;
            end
        end
        superCell = POP_STRUC.POPULATION(ind).superCell;
        fprintf(fp,'BestID: %4d;    SuperCEll:  %4d %4d %4d\n', ind, superCell);
        %%%%%%%%update Lattice%%%%%%%%%%%%%%
        lat_before = POP_STRUC.lat;
        if POP_STRUC.generation > 0
            [lat, lat_6] = UpdateH_tensor(lat_before, PressureTensor, ...
                ORG_STRUC.ExternalPressure, ORG_STRUC.GaussianWidth, ORG_STRUC.GaussianHeight,...
                POP_STRUC.lat_6);
        else
            [lat, lat_6] = UpdateH_tensor(lat_before, PressureTensor, ...
                ORG_STRUC.ExternalPressure, 0.5*ORG_STRUC.GaussianWidth, ORG_STRUC.GaussianHeight,...
                POP_STRUC.lat_6);
        end
        fprintf(fp,'\n');
        fprintf(fp,'    Lattice_before   ----->   Lattice_after            Delta H\n');
        fprintf(fp,'%6.3f  %6.3f  %6.3f;  %6.3f  %6.3f  %6.3f;  %6.3f  %6.3f  %6.3f\n',...
            lat_before(1,:),lat(1,:),lat(1,:)-lat_before(1,:));
        fprintf(fp,'%6.3f  %6.3f  %6.3f;  %6.3f  %6.3f  %6.3f;  %6.3f  %6.3f  %6.3f\n',...
            lat_before(2,:),lat(2,:),lat(2,:)-lat_before(2,:));
        fprintf(fp,'%6.3f  %6.3f  %6.3f;  %6.3f  %6.3f  %6.3f;  %6.3f  %6.3f  %6.3f\n',...
            lat_before(3,:),lat(3,:),lat(3,:)-lat_before(3,:));

        %%%%%%%%%End update Lattice%%%%%%%%%%%%%%5
        CreateCalcFolder();
        OFF_STRUC = POP_STRUC;
        OFF_STRUC.POPULATION = [];
        OFF_STRUC.lat   = lat;
        OFF_STRUC.lat_6 = lat_6;
        OFF_STRUC(1).POPULATION = POP_STRUC.POPULATION(1); %just make it has same structure
        OFF_STRUC.POPULATION(1) = POP_STRUC.POPULATION(ind);
        OFF_STRUC.POPULATION(1).howCome = 'keptBest';
        OFF_STRUC.POPULATION(1).Parents.parent = POP_STRUC.POPULATION(1).Number;
        OFF_STRUC.POPULATION(1).Parents.enthalpy = 0;
        OFF_STRUC.POPULATION(1).Step = [];

        next = 2; %The next ID number
        for i = 2 : ORG_STRUC.populationSize
            next = SoftMutation(next, ind, lat);
            if (next > ORG_STRUC.populationSize) || (next == 0)
                break; % exit if the maximum number is achieived or out of soft modes
            end
        end

        if ORG_STRUC.varcomp
            %%%%%%%%%%%%%%  if the 'basic' cell is chosen as parent alternative:  %%%%%%%%%
            %%%%%%%%%%%%%%         we want to relax this basic cell first         %%%%%%%%%
            if (mod(POP_STRUC.generation, ORG_STRUC.useBasicCell) == ORG_STRUC.useBasicCell-1) && (sum(superCell) > 3)
                OFF_STRUC.POPULATION(end+1) = POP_STRUC.POPULATION(ind);
                OFF_STRUC.POPULATION(end).COORDINATES = ORG_STRUC.bestBasicStructure;
                OFF_STRUC.POPULATION(end).Step = [];
                OFF_STRUC.POPULATION(end).numIons = ORG_STRUC.numIons;
                OFF_STRUC.POPULATION(end).Softmode_num = 0;
                OFF_STRUC.POPULATION(end).superCell = [1 1 1];
                OFF_STRUC.basicStructureNumber = length(OFF_STRUC.POPULATION); % remember the basic structure that we relax
                fprintf(fp, 'Adding the basic cell, to increase diversity---------------\n');
            elseif (mod(POP_STRUC.generation, ORG_STRUC.useBasicCell) == 0) ...
                && (sum(POP_STRUC.POPULATION(1).superCell) > 3) ...
                && (ind < POP_STRUC.basicStructureNumber)
                disp(' ');
                disp('Softmutating the basic cell, to increase diversity');
                fprintf(fp, 'Softmutating the basic cell, to increase diversity---------------\n');
                disp(' ');
                next = length(OFF_STRUC.POPULATION) + 1;
                for i = 1 : ORG_STRUC.populationSize
                    next = SoftMutation(next, POP_STRUC.basicStructureNumber, lat);
                    if (next > length(OFF_STRUC.POPULATION)+ORG_STRUC.populationSize) || (next==0)
                        break % exit if the maximum number is achieived or out of soft modes
                    end
                end
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%% THE NEW POPULATION
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        OFF_STRUC.generation = POP_STRUC.generation + 1;
        POP_STRUC = OFF_STRUC;
        WriteGenerationStart();
        %checkPreOptFings();
        META_Start_POP();
        fprintf(fp, [alignLine( sprintf('Generation %d', POP_STRUC.generation) ) '\n'] );
        fprintf(fp,'\n');

        fprintf(fp, [alignLine('Constant Lattice') '\n']);

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
