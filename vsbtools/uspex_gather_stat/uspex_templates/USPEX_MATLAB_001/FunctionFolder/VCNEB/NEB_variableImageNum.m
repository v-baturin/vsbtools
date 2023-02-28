function variableOFF_STRUC = NEB_variableImageNum()


% ---- modified
global ORG_STRUC
global POP_STRUC
global OFF_STRUC


%--------     System Parameters      ------%

OLD_numImages = ORG_STRUC.numImages;

weight      = ORG_STRUC.weight;
sumIons     = sum(ORG_STRUC.numIons);
dimension   = ORG_STRUC.dimension;
numDimension=3*(sumIons+dimension);
toleranceFing = ORG_STRUC.toleranceFing*10;

% Update the Energy tag

if POP_STRUC.POPULATION(1).preEnthalpy < ORG_STRUC.E_IniImage
    ORG_STRUC.E_IniImage   = POP_STRUC.POPULATION(1).preEnthalpy;
    ORG_STRUC.fp_IniImage  = calcIndFingerPrint( POP_STRUC.POPULATION(1) );
end
if POP_STRUC.POPULATION(OLD_numImages).preEnthalpy < ORG_STRUC.E_FinImage
    ORG_STRUC.E_FinImage   = POP_STRUC.POPULATION(OLD_numImages).preEnthalpy;
    ORG_STRUC.fp_FinImage  = calcIndFingerPrint( POP_STRUC.POPULATION(OLD_numImages) );
end

fp_first2nd = calcIndFingerPrint( POP_STRUC.POPULATION(2) );
fp_last2nd  = calcIndFingerPrint( POP_STRUC.POPULATION(OLD_numImages-1) );


dE_ini  = POP_STRUC.POPULATION(1).preEnthalpy             - ORG_STRUC.E_IniImage ;
dE_fin  = POP_STRUC.POPULATION(OLD_numImages).preEnthalpy - ORG_STRUC.E_FinImage;

fpD_first2nd =  cosineDistance( ORG_STRUC.fp_IniImage , fp_first2nd, weight );
fpD_last2nd  =  cosineDistance( ORG_STRUC.fp_FinImage , fp_last2nd,  weight );

%------------------------------------------%

dE = 0.002*sumIons;

maxLength     = ORG_STRUC.VarPathLength;
mapImageTable = zeros(OLD_numImages,2);

if ORG_STRUC.CalcType==1
    %-------------------------------%
    if (abs( dE_ini) < dE) && (fpD_first2nd<toleranceFing)
        startImage = 2;
        mapImageTable(1,:)=[2 1];
    else
        startImage = 1;
        mapImageTable(1,:)=[1 1];
    end
    if (abs( dE_fin ) < dE) && (fpD_last2nd<toleranceFing)
        endImage = OLD_numImages-1;
    else
        endImage = OLD_numImages;
    end
    %-------------------------------%
    k=2;
    for i = startImage+1:endImage
        %        disp('round : ')
        %        round( POP_STRUC.POPULATION(i).prePathLength/maxLength )
        if round( POP_STRUC.POPULATION(i).prePathLength/maxLength ) > 0
            mapImageTable(k,2) = mapImageTable(k-1,2) + round( POP_STRUC.POPULATION(i).prePathLength/maxLength ) ;
            mapImageTable(k,1) = i;
            k = k + 1;
        else
            mapImageTable(k-1,1) = i;
        end
    end
    
    lineTable=k-1;
    mapImageTable(1:k-1,:);
    NEW_numImages=mapImageTable(k-1,2);
else
    NEW_numImages=ORG_STRUC.numImages;
    mapImageTable(:,1)=1:OLD_numImages;
    mapImageTable(:,2)=1:OLD_numImages;
    lineTable=OLD_numImages;
end
%--------------------------------------------------------------------------


ORG_STRUC.numImages= NEW_numImages;
variableOFF_STRUC  = NEB_creatPOPStruct();

variableOFF_STRUC.genDone = 0;
variableOFF_STRUC.DoneOrder = zeros(NEW_numImages,1);
variableOFF_STRUC.bodyCount = 0;
variableOFF_STRUC.step      = OFF_STRUC.step;
variableOFF_STRUC.generation= OFF_STRUC.generation;
variableOFF_STRUC.resFolder = OFF_STRUC.resFolder;


for  i = 1:lineTable
    oldN = mapImageTable(i,1);
    newN = mapImageTable(i,2);
    
    variableOFF_STRUC.POPULATION(newN)=OFF_STRUC.POPULATION(oldN);
    if (oldN==2) && (newN==1)
        variableOFF_STRUC.POPULATION(newN).Fire.Ncount= 0;
        variableOFF_STRUC.POPULATION(newN).Fire.dt = ORG_STRUC.dt;
        variableOFF_STRUC.POPULATION(newN).Fire.a = 0.1;
        variableOFF_STRUC.POPULATION(newN).Fire.v = zeros(numDimension,1);
        variableOFF_STRUC.POPULATION(newN).Fire.P = 0;
    end
    if (endImage == OLD_numImages-1) && (oldN==OLD_numImages-1)
        variableOFF_STRUC.POPULATION(newN).Fire.Ncount= 0;
        variableOFF_STRUC.POPULATION(newN).Fire.dt = ORG_STRUC.dt;
        variableOFF_STRUC.POPULATION(newN).Fire.a = 0.1;
        variableOFF_STRUC.POPULATION(newN).Fire.v = zeros(numDimension,1);
        variableOFF_STRUC.POPULATION(newN).Fire.P = 0;
    end
    
end


%%%-----------mapImageTable
cd(ORG_STRUC.homePath);
disp(['Applying variable-Image-Method at VCNEB Step: ' num2str(OFF_STRUC.step)]);

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a');

if   0==mod(OFF_STRUC.step-1, ORG_STRUC.PrintStep) | (OFF_STRUC.step==ORG_STRUC.numSteps)
    fprintf(fp, '\n\n');
    
    fprintf(fp,'dE between %3d(first) and best initial Image : %10.4f (eV) \n',  1,              dE_ini);
    fprintf(fp,'dE between %3d(final)  and best final  Image : %10.4f (eV) \n',  OLD_numImages,  dE_fin );
    fprintf(fp,'dE tolerance : %10.3f (eV)\n', dE);
    fprintf(fp,'\n');
    fprintf(fp,'Finger Print Distance between %3d and best initial Image : %6.5f \n',  2,               fpD_first2nd );
    fprintf(fp,'Finger Print Distance between %3d and best final   Image : %6.5f \n',  OLD_numImages-1, fpD_last2nd  );
    fprintf(fp,'Finger Print tolerance : %6.3f \n', toleranceFing);
    
    if 1==ORG_STRUC.CalcType && ORG_STRUC.optVarImage==1
        fprintf(fp, '\nVariable Image Technic Mapping Table :\n');
        fprintf(fp, '\t\t  Images Sequential number at Step %3d ---> Step %3d \n',OFF_STRUC.step-1,OFF_STRUC.step );
        for iC = 1:size(mapImageTable,1)
            if mapImageTable(iC,1)>0
                fprintf(fp, '\t\t\t                     %3d     ->    %3d \n', mapImageTable(iC,:));
            end
        end
    end
end
for i = 1:OLD_numImages
    H(i)=OFF_STRUC.POPULATION(i).preEnthalpy;
end
fprintf(fp,'VCNEB Calculation Step %3d. \n', OFF_STRUC.step-1);
fprintf(fp, '\n   Barrier = %10.6f/%f eV/cell (->)/(<-).',   ...
    (max(H(2:OLD_numImages-1))-H(1)), (max(H(2:OLD_numImages-1))-H(OLD_numImages))  );
fprintf(fp, 'See file BarrierFig for barrier profile.\n');
fclose(fp);
%
fpath = [ORG_STRUC.resFolder, '/AuxiliaryFiles/ImageMapping.dat'];
if OFF_STRUC.step==1
    fp = fopen(fpath, 'w');
else
    fp = fopen(fpath, 'a');
end
fprintf(fp,'STEP %4d \n', OFF_STRUC.step-1);
for iC = 1:size(mapImageTable,1)
    if mapImageTable(iC,1)>0
        fprintf(fp, '%3d   %3d \n', mapImageTable(iC,:));
    end
end
fclose(fp);

%----------------------------------------------

numImages = NEW_numImages;

allImageVolume = zeros(numImages,1);
for i = 1:numImages
    allImageVolume(i)= det( variableOFF_STRUC.POPULATION(i).LATTICE );
end
for i = 2:numImages
    if allImageVolume(i) == 0
        prevImage=i-1; nextImage=0;
        for j=i+1:1:numImages
            if allImageVolume(j)>0
                nextImage=j;
                break;
            end
        end
        diffN = nextImage - prevImage;
        D_CARTECOORDS=variableOFF_STRUC.POPULATION(nextImage).CARTECOORDS-variableOFF_STRUC.POPULATION(prevImage).CARTECOORDS;
        D_LATTICE    =variableOFF_STRUC.POPULATION(nextImage).LATTICE  - variableOFF_STRUC.POPULATION(prevImage).LATTICE;
        
        variableOFF_STRUC.POPULATION(i).CARTECOORDS=variableOFF_STRUC.POPULATION(prevImage).CARTECOORDS+ D_CARTECOORDS/diffN;
        variableOFF_STRUC.POPULATION(i).LATTICE    =variableOFF_STRUC.POPULATION(prevImage).LATTICE   +    D_LATTICE/diffN;
        
    end
end


for i=1:numImages
    variableOFF_STRUC.POPULATION(i).COORDINATES = (  ( (variableOFF_STRUC.POPULATION(i).LATTICE')^(-1) )*variableOFF_STRUC.POPULATION(i).CARTECOORDS' )';
end

%----------------------

%
%
%

function FINGERPRINT = calcIndFingerPrint(POP)

global ORG_STRUC

LATTICE = POP.LATTICE;
numIons = POP.numIons;
COORDINATES = POP.COORDINATES;
atomType = ORG_STRUC.atomType;

[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, atomType);
[order, FINGERPRINT, atom_fing   ] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);

