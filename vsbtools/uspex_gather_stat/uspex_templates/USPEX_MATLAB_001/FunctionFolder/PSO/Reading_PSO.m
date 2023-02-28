function Reading_PSO(code, Ind_No, indic)

global POP_STRUC
global ORG_STRUC

cd (['CalcFold' num2str(indic)])

Gen = POP_STRUC.generation;
Step = POP_STRUC.POPULATION(Ind_No).Step;
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
maxErrors = ORG_STRUC.maxErrors;

Const_Lat = ORG_STRUC.constLattice; %fixed lattice, we don't optimize the cell
molecule  = ORG_STRUC.molecule; %For molecules, we need to reset the Molecule item (MOLCOORS, ZMATRIX, etc)
minDistMatrice = ORG_STRUC.minDistMatrice;

ID = ['Gen' num2str(Gen) '-Ind' num2str(Ind_No) '-Step' num2str(Step)]; %For ERROR output
TotalStep = length([ORG_STRUC.abinitioCode]);

% Check if there is error
GoodBad = Read_AbinitCode(code, 0, ID, Step);
if GoodBad==1
    [COORDINATES, LATTICE] = Read_Structure(code, Const_Lat);
    if ~distanceCheck(COORDINATES, LATTICE, numIons, minDistMatrice)
        POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
        POP_STRUC.POPULATION(Ind_No).Enthalpies(end) = 100000;
        %disp('The structure after relaxation cannot satisfy distance constraint');
        USPEXmessage(550, '', 0);
    end
    
    if POP_STRUC.POPULATION(Ind_No).Error <= maxErrors;
        POP_STRUC.POPULATION(Ind_No).COORDINATES = COORDINATES;
        POP_STRUC.POPULATION(Ind_No).LATTICE = LATTICE;
        POP_STRUC.POPULATION(Ind_No).Enthalpies(Step) = Read_AbinitCode(code, 1, ID, Step);
        if ((ORG_STRUC.optType == 6) | (ORG_STRUC.optType == 8)) & (Step == TotalStep)
            POP_STRUC.POPULATION(Ind_No).dielectric_tensor = Read_AbinitCode(code, 3, ID, Step);
        end
        
        if ((ORG_STRUC.optType == 7) & (Step == TotalStep)) | ((ORG_STRUC.optType == 8) & (Step == TotalStep-1))
            POP_STRUC.POPULATION(Ind_No).gap = Read_AbinitCode(code, 4, ID, Step);
        end
        
        if ORG_STRUC.optType == 9
            POP_STRUC.POPULATION(Ind_No).mag_moment = Read_AbinitCode(code, 5, ID, Step);
        end
        
        if (ORG_STRUC.optType == 8) & (Step == TotalStep-1)
            if (POP_STRUC.POPULATION(Ind_No).gap < 0.1) % got a metal here
                POP_STRUC.POPULATION(Ind_No).Step = Step + 2; % jump over the last step
            end
        end
        
        if ((ORG_STRUC.optType > 1100) & (ORG_STRUC.optType<1110)) & (Step == TotalStep)
            POP_STRUC.POPULATION(Ind_No).elasticMatrix=Read_AbinitCode(code, 6, ID, Step);
            elasticProperties=calcElasticProperties(POP_STRUC.POPULATION(Ind_No).elasticMatrix, ...
                POP_STRUC.POPULATION(Ind_No).numIons, ORG_STRUC.atomType, det(POP_STRUC.POPULATION(Ind_No).LATTICE));
            POP_STRUC.POPULATION(Ind_No).elasticProperties = elasticProperties;
        end
        
        if (ORG_STRUC.checkConnectivity == 1) & (Step == TotalStep)
            numIons = POP_STRUC.POPULATION(Ind_No).numIons;
            POP_STRUC.POPULATION(Ind_No).hardness = calcHardness(COORDINATES, LATTICE, numIons);
        end
        
        POP_STRUC.POPULATION(Ind_No).Step = Step + 1;
    end
    
    Clean_AbinitCode(code);
else
    if (code~=1) & (code~=8) % If error appears in the empirical code, we just skip this structure.
        POP_STRUC.POPULATION(Ind_No).Error = maxErrors + 1;
    else
        POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 1; %FOR VASP and QE, let's try again
    end
    unixCmd(['echo PROBLEM_reading Structure' num2str(Ind_No)]);
    unixCmd('pwd');
end
cd ..

