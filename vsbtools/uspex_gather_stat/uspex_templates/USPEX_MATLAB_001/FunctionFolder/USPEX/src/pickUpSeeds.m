function numOfGoods = pickUpSeeds()

% USPEX 10 (initially created by Guangrui Qian)
% Add seeds structure to the end of the POPULATION
% It has the following features

% 1: VASP-5.2 format
% 2: block varcomp support 
% 3: 

global POP_STRUC
global ORG_STRUC
warning off

%-------------------------------------
atomType   = ORG_STRUC.atomType;
numIons    = ORG_STRUC.numIons;
varcomp    = ORG_STRUC.varcomp;
molecule   = ORG_STRUC.molecule;
resFolder  = ORG_STRUC.resFolder;
STDMOL     = ORG_STRUC.STDMOL;
generation = POP_STRUC.generation;

cd Seeds

disp(' ');
disp('Read Seeds ... ')
disp(' ');

numOfGoods = 0;
% reading POSCARS
if exist([ 'POSCARS_', num2str(generation) ])
    if exist('POSCARS')
        unixCmd([ 'cat POSCARS >> POSCARS_' num2str(generation) ]);
    end
    unixCmd([ 'cp POSCARS_' num2str(generation) ' POSCARS' ]);
end

if exist('POSCARS')
    seedContent = {'THis is for the seed file'};
    [fid,message] = fopen('POSCARS');
    
    while 1
        tmp = fgetl(fid);      % system description
        if tmp==-1; break; end
        seedContent{end + 1} = tmp;
    end
else
    cd ..
    return;
end

%==========================================%

POSCARs     = splitSeedsContext( seedContent );

goodPOSCARs = postprocessSeeds( POSCARs, atomType, numIons, STDMOL);

disp(' ')
disp(' ')

numOfGoods = length(goodPOSCARs);


for ID = 1:length(goodPOSCARs)
    if generation == 1
        ORG_STRUC.initialPopSize = ORG_STRUC.initialPopSize + 1;
        Add = ORG_STRUC.initialPopSize;
    else
        Add = length(POP_STRUC.POPULATION) + 1;
    end
    
    POP_STRUC.POPULATION(Add).LATTICE     = goodPOSCARs(ID).lattice;
    POP_STRUC.POPULATION(Add).COORDINATES = goodPOSCARs(ID).coord;
    POP_STRUC.POPULATION(Add).numIons     = goodPOSCARs(ID).numIons;
    POP_STRUC.POPULATION(Add).howCome     = '  Seeds';
    
    if (ORG_STRUC.dimension==0) && (ORG_STRUC.varcomp==1)
         POP_STRUC.POPULATION(Add).INIT_FP = ...
                    fp_calc_001(goodPOSCARs(ID).lattice, goodPOSCARs(ID).coord, goodPOSCARs(ID).numIons);
    end
    
    if molecule==1
        POP_STRUC.POPULATION(Add).numMols = goodPOSCARs(ID).numMols;
        [type, MtypeLIST, numIons]=GetPOP_MOL(goodPOSCARs(ID).numMols);
        POP_STRUC.POPULATION(Add).typesAList = type;
        POP_STRUC.POPULATION(Add).MtypeLIST = MtypeLIST;
        readMOL(Add, 0);
    end
    if varcomp==1
        POP_STRUC.POPULATION(Add).numBlocks = goodPOSCARs(ID).numBlocks;
    end
    disp(['seed number ' num2str(ID) ' has been successfully added']);
end

disp(' ')
disp('End of pickup Seeds');
SeedsFile = ['../' resFolder '/Seeds_history'];
unixCmd([ 'echo Generation:' num2str(generation) '>> ' SeedsFile ]);
unixCmd([ 'cat POSCARS >> ../' SeedsFile ]);
unixCmd([ 'mv POSCARS POSCARS_' num2str(generation) ]);

cd ..
