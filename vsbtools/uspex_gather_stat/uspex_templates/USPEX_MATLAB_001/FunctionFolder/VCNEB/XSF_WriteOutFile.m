function XSF_WriteOutFile(step)


%--------------------------------------------------------------------%
global POP_STRUC
global ORG_STRUC

%--------------------------------------------------------------------%
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);
numIons  = ORG_STRUC.numIons;

%LATTICE
para = zeros(3,3);
pos = [ 0, 0, 0];


atomicNumber = zeros(1,sumIons);

iC=1;
for typeA = 1 : length(numIons)
    for i = 1 : numIons(typeA)
        atomicNumber(iC)=ORG_STRUC.atomType(typeA);
        iC = iC + 1;
    end
end

%-----------------------------------------------------
cd(ORG_STRUC.homePath);
fpath = [ORG_STRUC.resFolder '/PATH/path' num2str(POP_STRUC.step) '.xsf'];
fp = fopen(fpath, 'w');


%-- ASXF HEAD Format
fprintf(fp,'ANIMSTEPS %3d\n', numImages);
fprintf(fp,'CRYSTAL\n');

for i = 1:numImages
    fprintf(fp, ' PRIMVEC %2d\n', i);
    para = POP_STRUC.POPULATION(i).LATTICE;
    for j = 1:3
        fprintf(fp, '    %10.6f   %10.6f  %10.6f\n', para(j,:));
    end
    fprintf(fp, ' PRIMCOORD %2d\n', i );
    fprintf(fp, ' %2d %2d\n', sumIons, 1);
    for j = 1:sumIons
        fprintf(fp, ' %2d   %10.6f   %10.6f  %10.6f\n', atomicNumber(j), POP_STRUC.POPULATION(i).CARTECOORDS(j,:) );
    end
end

fclose(fp);
