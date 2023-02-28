function NEB_updateUSPEXFile()

global POP_STRUC
global ORG_STRUC


%
% Calculate the FingerPrint for all POP_STRUC
%
calcIndFingerPrint();

USPEXFile = [ORG_STRUC.homePath,'/' ORG_STRUC.resFolder '/USPEX.mat'];
if ~exist(USPEXFile,'file')
    %
    % Create the USPEX_STRUC
    %
    USPEX_STRUC.POPULATION = POP_STRUC.POPULATION;
    USPEX_STRUC.endPhases(1) = POP_STRUC.POPULATION(1);
    USPEX_STRUC.endPhases(2) = POP_STRUC.POPULATION(end);
    %
    %
    USPEX_STRUC.weight = ORG_STRUC.weight;
    USPEX_STRUC.PES    = calcNEBPES( USPEX_STRUC, 1 );
else
    load(USPEXFile);
    lenUSPEX = length( USPEX_STRUC.POPULATION );
    lenPOP   = length( POP_STRUC.POPULATION );
    %
    % Update the USPEX_STRUC
    %
    USPEX_STRUC.POPULATION(lenUSPEX+1:lenUSPEX+lenPOP) = POP_STRUC.POPULATION;
    USPEX_STRUC.endPhases(1) = POP_STRUC.POPULATION(1);
    USPEX_STRUC.endPhases(2) = POP_STRUC.POPULATION(end);
    
    
    USPEX_STRUC.weight = ORG_STRUC.weight;
    USPEX_STRUC.PES    = calcNEBPES( USPEX_STRUC, lenUSPEX+1 );
    
end

safesave(USPEXFile, USPEX_STRUC);

%
%
%
function calcIndFingerPrint()

global POP_STRUC
global ORG_STRUC

for whichInd = 1:length( POP_STRUC.POPULATION )
    LATTICE = POP_STRUC.POPULATION(whichInd).LATTICE;
    numIons = POP_STRUC.POPULATION(whichInd).numIons;
    COORDINATES = POP_STRUC.POPULATION(whichInd).COORDINATES;
    atomType = ORG_STRUC.atomType;
    
    [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
    [order, FINGERPRINT, atom_fing   ] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
    
    POP_STRUC.POPULATION(whichInd).order =  order;
    POP_STRUC.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
    POP_STRUC.POPULATION(whichInd).struc_entr = structureQuasiEntropy(numIons, atom_fing, numIons/sum(numIons));
    POP_STRUC.POPULATION(whichInd).S_order    = StructureOrder(FINGERPRINT, V, numIons, ORG_STRUC.deltaFing, ORG_STRUC.weight);
    
end
%
%
%
%
function PES = calcNEBPES( USPEX_STRUC, len0 )

weight = USPEX_STRUC.weight;
lenAll = length(USPEX_STRUC.POPULATION);


USPEX_STRUC;
PES=0;
