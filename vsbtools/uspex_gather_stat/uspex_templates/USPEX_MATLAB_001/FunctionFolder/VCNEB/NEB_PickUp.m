function PickUp()

global ORG_STRUC
global POP_STRUC

if ORG_STRUC.pickUpFolder == 0
    ORG_STRUC.pickUpFolder = str2num(ORG_STRUC.resFolder(8:end));
end


cd (['results' num2str( ORG_STRUC.pickUpFolder)]);
cd (['STEP/step' num2str(ORG_STRUC.pickUpGen)]);
disp(pwd);
load('POP_STRUC.mat')

cd( ORG_STRUC.homePath )

ORG_STRUC.pickedUP = 1;

ORG_STRUC.resFolder
unixCmd(['echo "This VCNEB Calculations has been picked up from an older calculation" >>' ORG_STRUC.resFolder '/PICKED_UP!!']);
save ([ORG_STRUC.resFolder '/PickedUP_POP_STRUC.mat'],'POP_STRUC');

ORG_STRUC.numImages=length(POP_STRUC.POPULATION);
ORG_STRUC.numSteps = ORG_STRUC.numSteps + POP_STRUC.step;


NEB_writeOutput(-1);
