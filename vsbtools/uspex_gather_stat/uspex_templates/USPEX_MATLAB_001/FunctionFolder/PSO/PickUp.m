function PickUp()

global ORG_STRUC
global POP_STRUC
global POOL_STRUC

if ORG_STRUC.pickUpFolder == 0
    ORG_STRUC.pickUpFolder = str2num(ORG_STRUC.resFolder(8:end))-1;
end

cd (['results' num2str(ORG_STRUC.pickUpFolder)])
try
    load('ANTISEEDS.mat')
catch
end

if ORG_STRUC.pickUpGen == 0
    ORG_STRUC.pickUpGen = 1;
    while exist (['generation' num2str( ORG_STRUC.pickUpGen)]) == 7
        ORG_STRUC.pickUpGen = ORG_STRUC.pickUpGen + 1;
    end
    ORG_STRUC.pickUpGen = ORG_STRUC.pickUpGen - 1;
end

cd (['generation' num2str(ORG_STRUC.pickUpGen)])
unixCmd('pwd');
load('POP_STRUC.mat')
if exist('POOL_STRUC.mat')
    load('POOL_STRUC.mat')
end
ORG_STRUC.initialPopSize = length(POP_STRUC.POPULATION);

cd ../..

if ORG_STRUC.varcomp == 1
    N_T = size(ORG_STRUC.numIons,1);
    splitting = zeros(1,N_T);
    findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
end

POP_STRUC.resFolder = ORG_STRUC.resFolder;
ORG_STRUC.pickedUP = 1;
disp('This Calculations has been picked up from an older calculation');
save ([ORG_STRUC.resFolder '/PickedUP_POP_STRUC.mat'],'POP_STRUC')

ORG_STRUC.numGenerations = ORG_STRUC.numGenerations + POP_STRUC.generation;

