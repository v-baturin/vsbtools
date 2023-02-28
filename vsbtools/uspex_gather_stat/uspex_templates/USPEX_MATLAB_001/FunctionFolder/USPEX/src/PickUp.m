function PickUp()

global ORG_STRUC
global POP_STRUC
global POOL_STRUC
global USPEX_STRUC

% Create the composition file for 300:
createCompostion_SingleBlock();

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
disp(pwd);
load('POP_STRUC.mat')
load('USPEX.mat')
if exist('POOL.mat')
    load('POOL.mat')
end


if ORG_STRUC.fixRndSeed > 0
    rng( ORG_STRUC.fixRndSeed+POP_STRUC.generation, 'twister' );
end

cd ../..


ORG_STRUC.initialPopSize = length(POP_STRUC.POPULATION);
ORG_STRUC.pickUpNCount = POP_STRUC.bodyCount;

if ORG_STRUC.varcomp == 1 || ~isempty(ORG_STRUC.maxAt)
    N_T = size(ORG_STRUC.numIons,1);
    splitting = zeros(1,N_T);
    findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
end

POP_STRUC.resFolder = ORG_STRUC.resFolder;
ORG_STRUC.pickedUP = 1;
disp('This Calculations has been picked up from an older calculation');
safesave ([ORG_STRUC.resFolder '/PickedUP_POP_STRUC.mat'],POP_STRUC)

ORG_STRUC.numGenerations = ORG_STRUC.numGenerations + POP_STRUC.generation;

safesave([ORG_STRUC.resFolder '/USPEX.mat'], USPEX_STRUC);
try
    safesave([ORG_STRUC.resFolder '/POOL.mat'], POOL_STRUC);
catch
end


%% This function creates a missing variable firstGeneSplit after PickUp:
function createCompostion_SingleBlock()

global ORG_STRUC

if isempty(ORG_STRUC.minAt) || isempty(ORG_STRUC.maxAt)
    ORG_STRUC.minAt = sum(ORG_STRUC.numIons);
    ORG_STRUC.maxAt = sum(ORG_STRUC.numIons);
end

N_T = size(ORG_STRUC.numIons,1);
splitting = zeros(1,N_T);
findSplit_VC(N_T, 0, ORG_STRUC.minAt, ORG_STRUC.maxAt, splitting);
IPS = ORG_STRUC.initialPopSize;

fp = fopen('Seeds/compositions', 'w');
if exist('Seeds/Anti-compositions')
    unixCmd('mv Seeds/Anti-compositions Seeds/Anti-compositions-back');
end
for i=1:size(ORG_STRUC.firstGeneSplit,1)
    for j=1:size(ORG_STRUC.firstGeneSplit,2)
        fprintf(fp, '%4d', ORG_STRUC.firstGeneSplit(i,j));
    end
    fprintf(fp, '\n');
end
fclose(fp);
