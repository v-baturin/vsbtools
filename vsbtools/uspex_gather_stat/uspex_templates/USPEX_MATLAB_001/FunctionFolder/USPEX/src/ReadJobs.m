function ID_next = ReadJobs(ID_now, ID_max, indic)

global ORG_STRUC
global POP_STRUC

dimension  = ORG_STRUC.dimension;
molecule   = ORG_STRUC.molecule;
varcomp    = ORG_STRUC.varcomp;

whichInd = ID_now;
ID_next  = ID_now;
Step     = POP_STRUC.POPULATION(whichInd).Step;
%N_Step   = ORG_STRUC.conv_till;
N_Step   = length([ORG_STRUC.abinitioCode]);

calcType   = [num2str(dimension) '' num2str(molecule) '' num2str(varcomp)];
if str2num(calcType) < 0
    tmp = -1*str2num(calcType);
    calcType = ['M' num2str(tmp)];
end


if POP_STRUC.POPULATION(whichInd).JobID > 0
    disp(['Structure' num2str(whichInd) ' step' num2str(Step) ' at CalcFold' num2str(indic) ]);
    if (ORG_STRUC.platform > 0) || (ORG_STRUC.numParallelCalcs > 1)
       disp(['JobID=' num2str(POP_STRUC.POPULATION(whichInd).JobID) ]);
    end
    doneOr = checkStatusC(whichInd);
    if doneOr
       if (POP_STRUC.POPULATION(whichInd).JobID == 0.01) ...
        | (POP_STRUC.POPULATION(whichInd).JobID == 0.02) % survived structure
           POP_STRUC.POPULATION(whichInd).Step = N_Step + 1;
       else
           Error = Reading(ORG_STRUC.abinitioCode(Step),whichInd, indic);
           if (Error == 0) && (dimension < 3) && (dimension ~= 1)
               eval(['Cell_transformation_' calcType '(whichInd)']);
           end
       end

%        if ORG_STRUC.conv_till < length([ORG_STRUC.abinitioCode])
%        POP_STRUC.finalOptimization == 1 % property calculation
           maxStep=length(ORG_STRUC.abinitioCode);
%        else
%           maxStep=ORG_STRUC.conv_till;
%        end

       POP_STRUC.POPULATION(whichInd).JobID = 0;

       if POP_STRUC.POPULATION(whichInd).Error > ORG_STRUC.maxErrors
          POP_STRUC.POPULATION(whichInd).Done = 1;
          POP_STRUC.POPULATION(whichInd).ToDo = 0;
          POP_STRUC.POPULATION(whichInd).Folder=0;
          ID_next = ID_max + 1;
          POP_STRUC.CalcFold_max    = ID_next;
          POP_STRUC.CalcFold(indic) = ID_next;

       %symmetrized structure for property calculation
       elseif POP_STRUC.POPULATION(whichInd).Step == ORG_STRUC.conv_till+1 && ORG_STRUC.conv_till < length([ORG_STRUC.abinitioCode]) && ORG_STRUC.abinitioCode(ORG_STRUC.conv_till+1) ==1
          lattice = POP_STRUC.POPULATION(whichInd).LATTICE;
          coor    = POP_STRUC.POPULATION(whichInd).COORDINATES;
          numIons = POP_STRUC.POPULATION(whichInd).numIons;
          order   = POP_STRUC.POPULATION(whichInd).order;
          atomType  = ORG_STRUC.atomType;

          cd([ORG_STRUC.homePath '/CalcFoldTemp']);
          POP_STRUC.POPULATION(whichInd).symg = anasym_stokes(lattice, coor, numIons, atomType, ORG_STRUC.SGtolerance);
    %% copy symmetrized.cif to CalcFold for property calculation
          currentCalcPath = [ORG_STRUC.homePath '/CalcFold' num2str(indic)];
          SymmetrizedFile = ([ORG_STRUC.homePath '/CalcFoldTemp/symmetrized.cif']);
          copyfile(SymmetrizedFile, currentCalcPath);
          cd(currentCalcPath);
          Natoms = Cif2poscar();
          totalNatoms = sum(numIons);
          if Natoms == totalNatoms
          elseif ORG_STRUC.conv_till+1 < length(ORG_STRUC.abinitioCode)
             copyfile('CONTCAR_old', 'outputPOSCAR');
          end
          delete('symmetrized.cif');
          cd(ORG_STRUC.homePath);
       %elseif POP_STRUC.POPULATION(whichInd).Step > length ([ORG_STRUC.abinitioCode])
       elseif POP_STRUC.POPULATION(whichInd).Step > maxStep
          
          if POP_STRUC.finalOptimization == 0  
            disp('Relaxation is done.')
            disp(' ')
            if ORG_STRUC.dimension == 0
                POP_STRUC.POPULATION(whichInd).COORDINATES = ...
                    moveCluster(POP_STRUC.POPULATION(whichInd).LATTICE,...
                    POP_STRUC.POPULATION(whichInd).COORDINATES);
            end                
            eval(['Fp_analysis_' calcType '(whichInd)']);
            POP_STRUC.bodyCount = POP_STRUC.bodyCount + 1;
            POP_STRUC.POPULATION(whichInd).Number = POP_STRUC.bodyCount;
            POP_STRUC.DoneOrder(whichInd) = POP_STRUC.bodyCount;
            eval(['WriteIndividualOutput_' calcType '(whichInd)']);
            POP_STRUC.POPULATION(whichInd).Folder=0;
            POP_STRUC.POPULATION(whichInd).Done = 1;
            POP_STRUC.POPULATION(whichInd).ToDo = 0;
            ID_next = ID_max + 1;
            POP_STRUC.CalcFold_max    = ID_next;
            POP_STRUC.CalcFold(indic) = ID_next;
         else
         %%%% waiting for more coding         
         end

       end
       safesave ('Current_POP.mat', POP_STRUC)
    else
        ID_next = 0;
    end
end
