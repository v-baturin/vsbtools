function NEB_writeFiles(step)

global ORG_STRUC
global POP_STRUC
%--------     System Parameters      ------%
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);

%--------------------------------------------------------------%
if 0==mod(POP_STRUC.step, ORG_STRUC.PrintStep)
    if     1==ORG_STRUC.FormatType
        XSF_WriteOutFile;
    elseif 2==ORG_STRUC.FormatType
        STM4_WriteOutFile(0);
    elseif 3==ORG_STRUC.FormatType
        XYZ_WriteOutFile;
    end
end

STM4_WriteOutFile(1);

%----------------------------------------------------------------%

NEB_findSymmetry();
NEB_writeOutput(1);


cd(ORG_STRUC.homePath);
%-------------------------------------------------------------------------%
fpath = [ POP_STRUC.resFolder '/SpaceGroup' ];
fp = fopen(fpath, 'w');
if 1==ORG_STRUC.CalcType
    fprintf(fp,'VCNEB Calcualtion Step %4d\n\n', POP_STRUC.step);
else
    fprintf(fp,'Relaxation Calcualtion Step %4d\n\n', POP_STRUC.step);
end
fprintf(fp,' Image       SpaceGroup\n');
for i = 1:numImages
    fprintf(fp, ' Image-%2d : %4d (%s)\n', i, POP_STRUC.POPULATION(i).spacegroupNumber, spaceGroups(POP_STRUC.POPULATION(i).spacegroupNumber) );
end
fclose(fp);

NEB_plotE();

if POP_STRUC.step >0
    fpath = [ POP_STRUC.resFolder,'/Force' ];
    fp = fopen(fpath, 'w');
    if 1==ORG_STRUC.CalcType
        fprintf(fp,'VCNEB Calcualtion Step %4d\n\n', POP_STRUC.step);
    else
        fprintf(fp,'Relaxation Calcualtion Step %4d\n\n', POP_STRUC.step);
    end
    fprintf(fp,'Image                  Force RootMeanSquare                                  Force Maximum  \n');
    fprintf(fp,'              Lattice                    Atomic        ------      Lattice                   Atomic  \n');
    fprintf(fp,'      [Total  Project  Elastic] [Total  Project  Elastic]  [Total  Project  Elastic] [Total  Project  Elastic] \n');
    for i = 1:numImages
        errcFnebRms = POP_STRUC.POPULATION(i).errcFnebRms;
        errcFproRms = POP_STRUC.POPULATION(i).errcFproRms;
        errcFelaRms = POP_STRUC.POPULATION(i).errcFelaRms;
        erraFnebRms = POP_STRUC.POPULATION(i).erraFnebRms;
        erraFproRms = POP_STRUC.POPULATION(i).erraFproRms;
        erraFelaRms = POP_STRUC.POPULATION(i).erraFelaRms;
        
        errcFnebMax = POP_STRUC.POPULATION(i).errcFnebMax;
        errcFproMax = POP_STRUC.POPULATION(i).errcFproMax;
        errcFelaMax = POP_STRUC.POPULATION(i).errcFelaMax;
        erraFnebMax = POP_STRUC.POPULATION(i).erraFnebMax;
        erraFproMax = POP_STRUC.POPULATION(i).erraFproMax;
        erraFelaMax = POP_STRUC.POPULATION(i).erraFelaMax;
        
        fprintf(fp,' %2d : [%7.4f %7.4f %7.4f] [%7.4f %7.4f %7.4f]--[%7.4f %7.4f %7.4f] [%7.4f %7.4f %7.4f]\n',i, errcFnebRms, errcFproRms, errcFelaRms, erraFnebRms, erraFproRms, erraFelaRms, errcFnebMax, errcFproMax, errcFelaMax, erraFnebMax, erraFproMax, erraFelaMax );
    end
    fclose(fp);
end


fpath = [POP_STRUC.resFolder '/ImageStructure'];
fp = fopen(fpath, 'w');
if 1==ORG_STRUC.CalcType
    fprintf(fp,'VCNEB Calcualtion Step %4d\n\n', POP_STRUC.step);
else
    fprintf(fp,'Relaxation Calcualtion Step %4d\n\n', POP_STRUC.step);
end
fprintf(fp, ' Image    a        b        c      alpha    beta   gamma  \n');
for i = 1:numImages
    Lattice = latConverter(POP_STRUC.POPULATION(i).LATTICE);
    if size(Lattice,1) == 1
        Lattice = Lattice';
    end
    Lattice(4:6) = Lattice(4:6)*180/pi;
    fprintf(fp, '  %2d  %7.4f  %7.4f  %7.4f,  %7.3f %7.3f %7.3f\n', i, Lattice(1:3), Lattice(4:6));
end

fclose(fp);


%------------------------------------------------------------------------------------
fpath = [ORG_STRUC.resFolder, '/AuxiliaryFiles/enthalpy_all.dat'];
fpath1= [ORG_STRUC.resFolder, '/AuxiliaryFiles/enthalpies_nospace.dat'];
if POP_STRUC.step==0
    %    unixCmd(['mkdir -p ' ORG_STRUC.resFolder '/AuxiliaryFiles']);
    unixCmd(['rm ' ORG_STRUC.resFolder '/AuxiliaryFiles/ImageMapping.dat']);
    fp = fopen(fpath,  'w');
    fp1= fopen(fpath1, 'w');
else
    fp = fopen(fpath, 'a');
    fp1= fopen(fpath1, 'a');
end
fprintf(fp,'%4d  %2d ', POP_STRUC.step, numImages);
for i = 1:numImages
    fprintf(fp, '%12.6f ', POP_STRUC.POPULATION(i).Enthalpy);
    fprintf(fp1,'%12.6f\n', POP_STRUC.POPULATION(i).Enthalpy);
end
fprintf(fp, '\n');
fclose(fp);
fclose(fp1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STM4_WriteOutFile(2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fpath = [ORG_STRUC.resFolder, '/AuxiliaryFiles/ImageStructure.dat'];

if POP_STRUC.step==0
    fp = fopen(fpath, 'w');
else
    fp = fopen(fpath, 'a');
end
fprintf(fp,'STEP  %4d \n', POP_STRUC.step);
for i = 1:numImages
    Lattice = latConverter(POP_STRUC.POPULATION(i).LATTICE);
    if size(Lattice,1) == 1
        Lattice = Lattice';
    end
    Lattice(4:6) = Lattice(4:6)*180/pi;
    fprintf(fp, '  %2d  %7.4f  %7.4f  %7.4f %7.3f %7.3f %7.3f\n', i, Lattice(1:3), Lattice(4:6));
end
fclose(fp);


makeFigures();
%---------------------------------------------------------------%
