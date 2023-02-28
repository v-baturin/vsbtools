function Error = Reading(code, Ind_No, indic)
% $Rev: 1493 $
% $Author: hyniu $
% $Date: 2016-03-03 16:09:19 -0500 (Thu, 03 Mar 2016) $

%This rountine is used to read some necessary tags from the output of Abintio Calculation;
%1, Check is the Calculation correctly done
%2, Read crystal structure (coordinate + lattice)
%3, Read target properties (energy/pressure tensor/dielectric constant, etc)
%4, clean output files in case they will be read at the following cycle
%Update: Qiang Zhu (2013/10/03)
%Add a new input argument to represent the path
%redefine the Fitness fucntion
%lastly updated by Qiang Zhu (2014/02/18)
%lastly updated by Haiyang Niu (2016/2/22)

global POP_STRUC
global ORG_STRUC

Error = 0;

cd (['CalcFold' num2str(indic)])
ParetoR = abs(ORG_STRUC.paretoRanking);
Gen     = POP_STRUC.generation;
Step    = POP_STRUC.POPULATION(Ind_No).Step;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;

try
    sumIons = sum(numIons);
catch
    % Fix for proteins interface (ticket #109).
    sumIons = 0;
end

maxErrors = ORG_STRUC.maxErrors;
Const_Lat = ORG_STRUC.constLattice; %fixed lattice, we don't optimize the cell
molecule  = ORG_STRUC.molecule;     %For molecules
atomType  = ORG_STRUC.atomType;
spin      = ORG_STRUC.spin;
dimension = ORG_STRUC.dimension;
minDist   = ORG_STRUC.minDistMatrice;
optType   = ORG_STRUC.optType;
goodBond  = ORG_STRUC.goodBonds;
N_val     = ORG_STRUC.NvalElectrons;
  val     = ORG_STRUC.valences;
collectForces  = ORG_STRUC.collectForces;
TotalStep = length(ORG_STRUC.abinitioCode);
%For ERROR output
ID = ['Gen' num2str(Gen) '-Ind' num2str(Ind_No) '-Step' num2str(Step)]; 

%Step 1; Check if the calculation is correctly done
GoodBad = Read_AbinitCode(code, 0, ID, Step);  
if GoodBad==1
    if (dimension ~=-4)
        %Use the last relaxation structure information for property calculation
        if Step == TotalStep && ORG_STRUC.conv_till < TotalStep
          COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
          LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;
        else
          [COORDINATES, LATTICE] = Read_Structure(code, Const_Lat);
        end
    else
        PROTEINS_STRUC = Read_Tinker_Structure(POP_STRUC.backbone_atoms);
    end
    if (molecule ~=1) && (dimension ~=-4) %check if it is correct
       
       if ORG_STRUC.conv_till == length([ORG_STRUC.abinitioCode])
        if ~distanceCheck(COORDINATES, LATTICE, numIons, minDist)
            POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
            POP_STRUC.POPULATION(Ind_No).Enthalpies(end) = 100000;
            disp('The structure after relaxation cannot satisfy distance constraint, please see the Warning ...');
            
            newPOP.numIons     = numIons;
            newPOP.LATTICE     = LATTICE;
            newPOP.COORDINATES = COORDINATES;
            brokenPOSCARStr = createPOSCARStr(POP_STRUC.POPULATION(Ind_No), ...
                                              newPOP, atomType);
            USPEXmessage(0, brokenPOSCARStr, 1);
            Error = 1;
        end
       end
    end
    
    if POP_STRUC.POPULATION(Ind_No).Error <= maxErrors;
        if (dimension ~=-4)
            POP_STRUC.POPULATION(Ind_No).COORDINATES = COORDINATES;
            POP_STRUC.POPULATION(Ind_No).LATTICE = LATTICE;
        else
            POP_STRUC.POPULATION(Ind_No).ANGLES      = PROTEINS_STRUC.angles;
            POP_STRUC.POPULATION(Ind_No).RESIDUES    = PROTEINS_STRUC.residues;
            POP_STRUC.POPULATION(Ind_No).SEC_STRUCT  = PROTEINS_STRUC.secondary_structure;
            POP_STRUC.POPULATION(Ind_No).COORDINATES = PROTEINS_STRUC.backbone_crd_norm;
            POP_STRUC.POPULATION(Ind_No).LATTICE     = PROTEINS_STRUC.lattice;
            POP_STRUC.POPULATION(Ind_No).numIons     = PROTEINS_STRUC.numIons;
            
            % weight is always 1 for proteins since we treat alpha-carbon atoms only:
            ORG_STRUC.weight   = 1;
            % This is for weight calculation in KeepBestStructures():
            ORG_STRUC.numIons  = PROTEINS_STRUC.numIons;
            ORG_STRUC.atomType = PROTEINS_STRUC.atomType;
        end
        POP_STRUC.POPULATION(Ind_No).Enthalpies(Step) = Read_AbinitCode(code, 1, ID, Step);
	%% Sometimes symmetrization change the numIons so the enthalpy comes from vasp is different, 
	%% In this case we use the enthalpy in the previous step, Haiyang Niu and ZX
        if Step > ORG_STRUC.conv_till
            POP_STRUC.POPULATION(Ind_No).Enthalpies(Step) = POP_STRUC.POPULATION(Ind_No).Enthalpies(ORG_STRUC.conv_till);
        end
	%%------------------------------------------------
        if molecule == 1
            readMOL(Ind_No, code); %for MOLECULES
        end
        
        if spin == 1
            atomicMagmom = Read_AbinitCode(code, -5, ID, Step);
            magmom_ions  = magTypeIdentify([POP_STRUC.POPULATION(Ind_No).magmom_ions(1), atomicMagmom(1:end,end)'], POP_STRUC.POPULATION(Ind_No).magmom_ini);
            
            POP_STRUC.POPULATION(Ind_No).magmom_ions(Step,:) = magmom_ions;
            POP_STRUC.POPULATION(Ind_No).mag_moment = sum(magmom_ions(2:end));
            if Step<length([ORG_STRUC.abinitioCode])
                POP_STRUC.POPULATION(Ind_No).magmom_ions(Step+1,:)=POP_STRUC.POPULATION(Ind_No).magmom_ions(Step,:);
            end
        end
        if collectForces == 1
            if ~isfield(POP_STRUC.POPULATION(Ind_No), 'RelaxStep')
                POP_STRUC.POPULATION(Ind_No).RelaxStep(1).Energy = [];
                POP_STRUC.POPULATION(Ind_No).RelaxStep(1).LATTICE= [];
                POP_STRUC.POPULATION(Ind_No).RelaxStep(1).COORD = [];
                POP_STRUC.POPULATION(Ind_No).RelaxStep(1).FORCE = [];
                POP_STRUC.POPULATION(Ind_No).RelaxStep(1).STRESS= [];
            end
            allLattice=  Read_AbinitCode(code, 9, ID, Step);
            allCoords =  Read_AbinitCode(code, 8, ID, Step);
            allForces =  Read_AbinitCode(code, 7, ID, Step);
            allStess  =  Read_AbinitCode(code,-2, ID, Step);
            allEnergy =  Read_AbinitCode(code,-1, ID, Step);

            relaxStep = min( size(allCoords, 1), size(allForces, 1) )/sum(numIons);
            for k = 1:relaxStep
                if isempty(allLattice)
                    tmpLattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
                elseif 3 == size(allLattice, 1)
                    tmpLattice = allLattice;
                else
                    tmpLattice = allLattice( (k-1)*3+1:k*3, : );
                end
                POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).Energy(k)                   = allEnergy(k);
                POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).LATTICE(1:3, 1:3, k)        = tmpLattice;
                POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).COORD  (1:sumIons, 1:3, k)  = allCoords( (k-1)*sumIons+1:k*sumIons, : );
                POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).FORCE  (1:sumIons, 1:3, k)  = allForces( (k-1)*sumIons+1:k*sumIons, : );
                POP_STRUC.POPULATION(Ind_No).RelaxStep(Step).STRESS (1:3, 1:3, k)        = allStess(  (k-1)*3+1:k*3, : );
            end
        end
        
        if ((optType == 6) || (optType == 11) || (optType == 8) ...
             || (ParetoR == 6) || (ParetoR == 11) || (ParetoR == 8)) && (Step == TotalStep)
            POP_STRUC.POPULATION(Ind_No).dielectric_tensor = Read_AbinitCode(code, 3, ID, Step);
            POP_STRUC.POPULATION(Ind_No).birefringence = Read_AbinitCode(code, 10, ID, Step);
        end
        
        if (((optType == 7) || (ParetoR == 7)) && (Step == TotalStep)) ...
           || (((optType == 8) || (ParetoR == 8)) && (Step == TotalStep-1))
            POP_STRUC.POPULATION(Ind_No).gap = Read_AbinitCode(code, 4, ID, Step);
        end
        
        if optType == 9 || ParetoR == 9
            POP_STRUC.POPULATION(Ind_No).mag_moment = Read_AbinitCode(code, 5, ID, Step);
        end
        
        if ((optType == 8) || (ParetoR == 8))&& (Step == TotalStep-1)
            if (POP_STRUC.POPULATION(Ind_No).gap < 0.1) % got a metal here
                POP_STRUC.POPULATION(Ind_No).Step = Step + 2; % Skip it
            end
        end

        if (optType == 14 || ParetoR == 14)  && (Step == TotalStep)
            % flag should be set according to TE_goal
            POP_STRUC.POPULATION(Ind_No).TE_property = Read_AbinitCode(code, 3, ID, Step);
        end
	if (optType ==17 || ParetoR == 17) && (Step == TotalStep)
	    POP_STRUC.POPULATION(Ind_No).Fphon = Read_AbinitCode(code, 3, ID, Step);
	end 

        if (((optType > 1100) && (optType<1112)) || ((ParetoR > 1100) && (ParetoR < 1112)))  && (Step == TotalStep) 
            elasticMatrix     = Read_AbinitCode(code, 6, ID, Step);
            elasticProperties = calcElasticProperties(elasticMatrix, ...
                                     numIons, atomType, det(LATTICE));
            POP_STRUC.POPULATION(Ind_No).elasticProperties = elasticProperties;
            POP_STRUC.POPULATION(Ind_No).elasticMatrix     = elasticMatrix;
        end
        
        if (ORG_STRUC.checkConnectivity == 1) && (Step == TotalStep)
            POP_STRUC.POPULATION(Ind_No).hardness = calcHardness(LATTICE,...
            COORDINATES, numIons, atomType, goodBond, N_val, val, dimension);
        end
        
        POP_STRUC.POPULATION(Ind_No).Step = Step + 1;

        % Updating step must be behind the above line
        % or else, the update to Step will be overwritten.
        if ((optType == 14) && (Step == TotalStep - 1))
          % check completeness of VASP results
          is_vasp_finished = Read_AbinitCode(code, -100, ID, Step);
          if 0 == is_vasp_finished
            % the structure is too bad, drop it.
            [nothing, results] = unix('wc -l EIGENVAL');
            disp(['Number of lines in EIGENVAL: ' results]);
            POP_STRUC.POPULATION(Ind_No).Step = Step + 2; % skip it
          end
        end

    end
    
    Clean_AbinitCode(code);
    
else
    Error = 1;
    if (code~=1) && (code~=8) && (code~=9) && (code~=14) && (code~=17) 
    % If error appears in the empirical code, skip this structure.
        POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
    else
    % For VASP/QE/FHI, let's try again 
        POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 1; 
    end
    [a,b]=unix(['echo PROBLEM_reading Structure' num2str(Ind_No)]);
    disp(pwd);
end
cd ..


%=====================
function brokenPOSCARStr = createPOSCARStr(oldPOP, POP, atomType )

atomTypStr = [];
for i = 1:length(atomType)
    atomTypStr = [atomTypStr, ' ', megaDoof( atomType(i) ) ' '];
end


% after Relaxation
brokenPOSCARStr= ['The structure after relaxation cannot satisfy distance constraint\n', '1.00\n'];

for i = 1:3
    brokenPOSCARStr= [ brokenPOSCARStr, num2str(POP.LATTICE(i,:)) '\n'];
end

brokenPOSCARStr= [ brokenPOSCARStr, atomTypStr, '\n'];
brokenPOSCARStr= [ brokenPOSCARStr, num2str(POP.numIons),'\n', 'Direct\n'];
for i = 1:size( POP.COORDINATES, 1 )
    brokenPOSCARStr= [ brokenPOSCARStr, num2str(POP.COORDINATES(i,:)) '\n' ];
end

% The original Structure
brokenPOSCARStr= [brokenPOSCARStr, '\n\nStruture before relaxation, for checking\n', '1.00\n'];

for i = 1:3
    brokenPOSCARStr= [ brokenPOSCARStr, num2str(oldPOP.LATTICE(i,:)) '\n'];
end

brokenPOSCARStr= [ brokenPOSCARStr, atomTypStr, '\n'];
brokenPOSCARStr= [ brokenPOSCARStr, num2str(oldPOP.numIons),'\n', 'Direct\n'];
for i = 1:size( oldPOP.COORDINATES, 1 )
    brokenPOSCARStr= [ brokenPOSCARStr, num2str(oldPOP.COORDINATES(i,:)) '\n' ];
end
