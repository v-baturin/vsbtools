function SubmitJobs_PSO()
global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    DO_NOW = 0;
    if isempty (whichInd)
        stillLeft = find([POP_STRUC.POPULATION(:).ToDo]);
        if ~isempty(stillLeft)
            DO_NOW = stillLeft(1);
            POP_STRUC.POPULATION(DO_NOW).Folder = indic;
            POP_STRUC.POPULATION(DO_NOW).ToDo = 0;
        end
    elseif POP_STRUC.POPULATION(whichInd).JobID == 0
        DO_NOW = whichInd;
    end
    if DO_NOW
        skipOld = 1;
        if ORG_STRUC.reoptOld == 1
            skipOld = 0;
        elseif isempty(POP_STRUC.POPULATION(DO_NOW).Parents)
            skipOld = 0;
        elseif ~(isfield(POP_STRUC.POPULATION(DO_NOW).Parents, 'keptAsBest') || isfield(POP_STRUC.POPULATION(DO_NOW).Parents, 'keptAsConvexHull'))
            skipOld = 0;
        end
        Step = POP_STRUC.POPULATION(DO_NOW).Step;
        
        %between all steps except last before last and the last step, all except kept as best and keptAsConvexHull
        perturbIt = 1;
        if ~(isempty(POP_STRUC.POPULATION(DO_NOW).Parents))
            if (isfield(POP_STRUC.POPULATION(DO_NOW).Parents, 'keptAsBest') || isfield(POP_STRUC.POPULATION(DO_NOW).Parents, 'keptAsConvexHull'))
                if ORG_STRUC.reoptOld == 0 % do not shake the survived structures, when reoptOld = 0
                    perturbIt = 0;
                end
            end
        end
        if (Step > 1) && (Step < length(ORG_STRUC.abinitioCode)-1)
            if perturbIt
                POP_STRUC.POPULATION(DO_NOW).COORDINATES = perturbCoords(POP_STRUC.POPULATION(DO_NOW).COORDINATES, POP_STRUC.POPULATION(DO_NOW).LATTICE,1);
            end
        end
        if skipOld
            POP_STRUC.POPULATION(DO_NOW).JobID = 0.01;
        elseif ORG_STRUC.abinitioCode(Step) == 0   % no optimization at all! (used in order optimization)
            POP_STRUC.POPULATION(DO_NOW).JobID = 0.02;
        else
            cd (['CalcFold' num2str(indic)])
            Write_AbinitCode(ORG_STRUC.abinitioCode(Step), DO_NOW);
            POP_STRUC.POPULATION(DO_NOW).JobID = submitJob(DO_NOW);
            cd ..
        end
        save ('Current_POP.mat', 'POP_STRUC')
    end
end
