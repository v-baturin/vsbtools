function ReadJobs_PSO()
global ORG_STRUC
global POP_STRUC

for indic = 1:ORG_STRUC.numParallelCalcs
    % check whether calculations are still performed and read out the output if they are done
    whichInd = find([POP_STRUC.POPULATION(:).Folder]==indic);
    if ~isempty (whichInd)
        Step = POP_STRUC.POPULATION(whichInd).Step;
        disp(['Structure' num2str(whichInd) ' step' num2str(Step) ' at CalcFold' num2str(indic) ]);
        if POP_STRUC.POPULATION(whichInd).JobID
            if ORG_STRUC.platform > 0
                disp(['JobID=' num2str(POP_STRUC.POPULATION(whichInd).JobID) ]);
            end
            doneOr = checkStatusC(whichInd);
            if doneOr
                if POP_STRUC.POPULATION(whichInd).JobID == 0.01   % reoptOld = 0, survived structure
                    POP_STRUC.POPULATION(whichInd).Step = length([ORG_STRUC.abinitioCode]) + 1;
                elseif POP_STRUC.POPULATION(whichInd).JobID ~= 0.02   % no optimization
                    Reading_PSO(ORG_STRUC.abinitioCode(Step),whichInd, indic);
                end
                
                POP_STRUC.POPULATION(whichInd).JobID = 0;
                
                if POP_STRUC.POPULATION(whichInd).Error > ORG_STRUC.maxErrors
                    POP_STRUC.POPULATION(whichInd).Done = 1;
                    POP_STRUC.POPULATION(whichInd).ToDo = 0;
                    POP_STRUC.POPULATION(whichInd).Folder=0;
                    
                elseif POP_STRUC.POPULATION(whichInd).Step > length ([ORG_STRUC.abinitioCode])
                    POP_STRUC.POPULATION(whichInd).Done = 1;
                    POP_STRUC.POPULATION(whichInd).ToDo = 0;
                    POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
                    POP_STRUC.POPULATION(whichInd).Number = POP_STRUC.bodyCount;
                    LATTICE = POP_STRUC.POPULATION(whichInd).LATTICE;
                    numIons = POP_STRUC.POPULATION(whichInd).numIons;
                    COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
                    atomType = ORG_STRUC.atomType;
                    [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
                    [order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix,typ_i, typ_j, numIons);
                    POP_STRUC.POPULATION(whichInd).order =  order;
                    POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
                    disp('Relaxation is done.')
                    disp(' ')
                    POP_STRUC.POPULATION(whichInd).struc_entr = structureQuasiEntropy(numIons, atom_fing, numIons/sum(numIons));
                    POP_STRUC.POPULATION(whichInd).S_order    = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);
                    
                    POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
                    WriteIndividualOutput_PSO(whichInd);
                    POP_STRUC.POPULATION(whichInd).Folder=0;
                end
                save ('Current_POP.mat', 'POP_STRUC')
            end
        end
    end
end
