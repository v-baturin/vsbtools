function NEB_generateNewImage()

% ---- modified
global ORG_STRUC
global POP_STRUC
global OFF_STRUC

%--------     System Parameters      ------%
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);
dimension  = ORG_STRUC.dimension;
numDimension=3*(sumIons+dimension);

%------------------------------------------------------------------------%
para = zeros(3);
para0= zeros(3);

%------------------------------------------------------------------------%

for i = 1:numImages
    if 0 == POP_STRUC.POPULATION(i).freezing
        tmp = ( reshape( POP_STRUC.POPULATION(i).X_move',3,sumIons+dimension ) )';  %-- need to check
        % u_move reshape to 3x(sumIons+dimension) Matrix
        %    display('reshape of X_move');
        %    tmp
        POP_STRUC.POPULATION(i).CARTECOORDS_move = tmp(1+dimension:sumIons+dimension, :);
        
        POP_STRUC.POPULATION(i).LATTICE_move    = tmp(1:dimension,:);
        
        %-----------  The 1st atom is fixed during the relax!!!-------------------------------------------------------
        for j = sumIons:-1:1
            POP_STRUC.POPULATION(i).CARTECOORDS_move(j,:) = POP_STRUC.POPULATION(i).CARTECOORDS_move(j,:) - POP_STRUC.POPULATION(i).CARTECOORDS_move(1,:);
        end
        L = POP_STRUC.POPULATION(i).LATTICE;
        POP_STRUC.POPULATION(i).EPSILON_move   = ( (L')^(-1)*(POP_STRUC.POPULATION(i).LATTICE_move) );
        POP_STRUC.POPULATION(i).COORDINATES_move= ( (L')^(-1)*(POP_STRUC.POPULATION(i).CARTECOORDS_move') )';
        
    else
        POP_STRUC.POPULATION(i).LATTICE_move = zeros(3,3);
        POP_STRUC.POPULATION(i).EPSILON_move = zeros(3,3) ;
        POP_STRUC.POPULATION(i).COORDINATES_move = zeros(sumIons,3) ;
    end
end



%%---------------------------------------------------
%disp(' New Image Cell & Atom Displacement')
%POP_STRUC.POPULATION(i).LATTICE_move
%POP_STRUC.POPULATION(i).EPSILON_move
%POP_STRUC.POPULATION(i).COORDINATES_move
%POP_STRUC.POPULATION(i).CARTECOORDS_move

%---------------------------------------------------


OFF_STRUC = NEB_creatPOPStruct();

OFF_STRUC.genDone = 0;
OFF_STRUC.DoneOrder = zeros(numImages,1);
OFF_STRUC.bodyCount = 0;
OFF_STRUC.step = POP_STRUC.step+1;
OFF_STRUC.generation = POP_STRUC.generation+1;
OFF_STRUC.resFolder = POP_STRUC.resFolder;

%NEB_checkMove();

for i=1:numImages
    L = POP_STRUC.POPULATION(i).LATTICE ;
    
    currNewCellLat = POP_STRUC.POPULATION(i).LATTICE_move + POP_STRUC.POPULATION(i).LATTICE;
    currNewCellEps = POP_STRUC.POPULATION(i).EPSILON_move + POP_STRUC.POPULATION(i).EPSILON;
    OFF_STRUC.POPULATION(i).LATTICE = round( (currNewCellLat)*1e+12 )/1e+12;
    OFF_STRUC.POPULATION(i).EPSILON = round( (currNewCellEps)*1e+12 )/1e+12;
    %-------------------------
    tmp = (POP_STRUC.POPULATION(i).CARTECOORDS_move + POP_STRUC.POPULATION(i).CARTECOORDS);
    OFF_STRUC.POPULATION(i).COORDINATES = ( round( (L')^(-1)*(tmp')*1e+12 )/1e+12 )' ;
    OFF_STRUC.POPULATION(i).CARTECOORDS = round( OFF_STRUC.POPULATION(i).COORDINATES*OFF_STRUC.POPULATION(i).LATTICE*1e+12 )/1e+12;
    OFF_STRUC.POPULATION(i).Fire = POP_STRUC.POPULATION(i).Fire;
    
    OFF_STRUC.POPULATION(i).VOLUME=abs( det( OFF_STRUC.POPULATION(i).LATTICE ) );
    OFF_STRUC.POPULATION(i).preEnthalpy  = POP_STRUC.POPULATION(i).Enthalpy;
    OFF_STRUC.POPULATION(i).prePathLength  = POP_STRUC.POPULATION(i).pathLength;
    OFF_STRUC.POPULATION(i).prePress  = POP_STRUC.POPULATION(i).Press;
    OFF_STRUC.POPULATION(i).preVOLUME = POP_STRUC.POPULATION(i).VOLUME;
    OFF_STRUC.POPULATION(i).Error=0;
    if POP_STRUC.POPULATION(i).freezing==1
        OFF_STRUC.POPULATION(i)=POP_STRUC.POPULATION(i);
        OFF_STRUC.POPULATION(i).JobID=-1;
        OFF_STRUC.POPULATION(i).ToDo = 0;
        OFF_STRUC.POPULATION(i).Done = 1;
    end
end

%NEB_bulkCalc();

NEB_displayCalculation();


%------------------------------------------------%
%ORG_STRUC.optVarImage=1;
if ORG_STRUC.optVarImage==0 || POP_STRUC.step==0 || ORG_STRUC.CalcType==2
    POP_STRUC = OFF_STRUC;
    NEB_pathLength();
else
    POP_STRUC = OFF_STRUC;
    NEB_pathLength();
    variableOFF_STRUC = NEB_variableImageNum();
    POP_STRUC = variableOFF_STRUC;
end

%disp(['in NEB_generateNEWImage.m: length of Image ' num2str(length(POP_STRUC.POPULATION))])

%save ('.CHECK_POP.mat', POP_STRUC);
%save ('.CHECK_OFF.mat', OFF_STRUC);
