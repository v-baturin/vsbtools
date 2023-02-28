function PickUp()

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

if ORG_STRUC.pickUpFolder == 0
    ORG_STRUC.pickUpFolder = str2num(ORG_STRUC.resFolder(8:end))-1;
end

cd (['results' num2str(ORG_STRUC.pickUpFolder)])

if ORG_STRUC.pickUpGen == 0
    ORG_STRUC.pickUpGen = 1;
    while exist (['generation' num2str( ORG_STRUC.pickUpGen)]) == 7
        ORG_STRUC.pickUpGen = ORG_STRUC.pickUpGen + 1;
    end
    ORG_STRUC.pickUpGen = ORG_STRUC.pickUpGen - 1;
end

cd (['generation' num2str(ORG_STRUC.pickUpGen)])
disp(pwd);
load('POP_STRUC.mat')
load('USPEX.mat')
cd ../..

ORG_STRUC.initialPopSize = length(POP_STRUC.POPULATION);

POP_STRUC.resFolder = ORG_STRUC.resFolder;
ORG_STRUC.pickedUP = 1;
disp('This Calculations has been picked up from an older calculation');
safesave ([ORG_STRUC.resFolder '/PickedUP_POP_STRUC.mat'],POP_STRUC)
ORG_STRUC.numGenerations = ORG_STRUC.numGenerations + POP_STRUC.generation;

safesave([ORG_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
