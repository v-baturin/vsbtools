function META_ReadJobs()
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
              META_Reading(ORG_STRUC.abinitioCode(Step),whichInd, indic);
              if POP_STRUC.finalOptimization == 1  % 'final' optimization 
                 maxStep = length([ORG_STRUC.abinitioCode]);
              else
                 maxStep = ORG_STRUC.conv_till;
              end

              POP_STRUC.POPULATION(whichInd).JobID = 0;

              if POP_STRUC.POPULATION(whichInd).Error > ORG_STRUC.maxErrors
                  POP_STRUC.POPULATION(whichInd).Done = 1;
                  POP_STRUC.POPULATION(whichInd).ToDo = 0;
                  POP_STRUC.POPULATION(whichInd).Folder=0;
              %Note here Step might be updated in Reading.m
              elseif POP_STRUC.POPULATION(whichInd).Step > maxStep  
                  numIons = POP_STRUC.POPULATION(whichInd).numIons;
                  LATTICE = POP_STRUC.POPULATION(whichInd).LATTICE;
                  COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
                  atomType = ORG_STRUC.atomType;
                  [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
                  [order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
                  POP_STRUC.POPULATION(whichInd).order =  order;
                  POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;

                  if POP_STRUC.finalOptimization == 0
                     POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
                     POP_STRUC.POPULATION(whichInd).Number = POP_STRUC.bodyCount;
                     POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
                     POP_STRUC.POPULATION(whichInd).lat0    = LATTICE;
                     POP_STRUC.POPULATION(whichInd).coords0 = COORDINATES;
                     POP_STRUC.POPULATION(whichInd).PressureTensor0 = POP_STRUC.POPULATION(whichInd).PressureTensor;
                     META_WriteIndividualOutput(whichInd, Step, 1);  %WriteIndividuals before full relaxation
                     disp('Constant lattice relaxation is done.')
                  else

                     META_WriteIndividualOutput(whichInd, Step, 2);  %WriteIndividuals after  full relaxation
                     disp('Full relaxation mode is done.')
                  end
                  disp(' ')
                  POP_STRUC.POPULATION(whichInd).Folder=0;
                  POP_STRUC.POPULATION(whichInd).Done = 1;
                  POP_STRUC.POPULATION(whichInd).ToDo = 0;
              end
              safesave ('Current_POP.mat', POP_STRUC)
           end
       end
    end
end
