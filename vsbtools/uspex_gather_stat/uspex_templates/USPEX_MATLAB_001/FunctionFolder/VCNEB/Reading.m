function Error = Reading(code, Ind_No, indic)


%This rountine is used to read some necessary tags from the output of Abintio Calculation;
%1, Check is the Calculation correctly done
%2, Read crystal structure (coordinate + lattice)
%3, Read target properties (energy/pressure tensor/dielectric constant, etc)
%4, clean output files in case they will be read at the following cycle

global POP_STRUC
global ORG_STRUC

Error = 0;

cd ([ ORG_STRUC.homePath '/CalcFold' num2str(indic)])

Gen = POP_STRUC.generation;
Step = POP_STRUC.POPULATION(Ind_No).Step;
numIons = ORG_STRUC.numIons;
maxErrors = 3;%ORG_STRUC.maxErrors;

Const_Lat = 1;%ORG_STRUC.constLattice; %fixed lattice, we don't optimize the cell

ID = ['Gen' num2str(Gen) '-Ind' num2str(Ind_No) '-Step' num2str(Step)]; %For ERROR output
TotalStep = length([ORG_STRUC.abinitioCode]);



GoodBad = Read_AbinitCode(code, 0, ID, Step);  %Step 1; Check if the calculation is correctly done
if GoodBad==1
    %    [COORDINATES, LATTICE] = Read_Structure(code, Const_Lat);
    
    if POP_STRUC.POPULATION(Ind_No).Error <= maxErrors;
        %        POP_STRUC.POPULATION(Ind_No).COORDINATES = COORDINATES;
        %        POP_STRUC.POPULATION(Ind_No).LATTICE = LATTICE;
        
        POP_STRUC.POPULATION(Ind_No).Enthalpy(Step)   = Read_AbinitCode(code, 1, ID, Step);
        POP_STRUC.POPULATION(Ind_No).cellStressMatirx = Read_AbinitCode(code, 2, ID, Step);
        origReadingForce    = Read_AbinitCode(code, 7, ID, Step);
        
        
        POP_STRUC.POPULATION(Ind_No).atomForcesMatrix = origReadingForce;
        if     code==1  %--> VASP
            
        elseif code==3  %--> GULP
            L=POP_STRUC.POPULATION(Ind_No).LATTICE;
            POP_STRUC.POPULATION(Ind_No).atomForcesMatrix = -( (L'*L)*L^(-1)*origReadingForce' )'; %eV/A
        elseif code==7  %-->CP2K
            
        elseif code==8  %-->QE
            
        end
        
        POP_STRUC.POPULATION(Ind_No).Step = Step + 1;
        POP_STRUC.POPULATION(Ind_No).Error=0;
    end
    
    Clean_AbinitCode(code);
    
else
    Error = 1;
    if (code~=1) && (code~=8) && (code~=9)  % If error appears in the empirical code, we just skip this structure.
        POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
    else
        POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 1; %FOR VASP/QE/FHI, let's try again
    end
    disp(['PROBLEM_reading Structure ' num2str(Ind_No)]);
    disp(pwd);
end
cd ..
