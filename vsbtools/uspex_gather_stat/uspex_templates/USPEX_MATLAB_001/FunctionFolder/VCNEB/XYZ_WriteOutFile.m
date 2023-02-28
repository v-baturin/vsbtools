function XYZ_WriteOutFile(step)


%--------------------------------------------------------------------%
global POP_STRUC
global ORG_STRUC

%--------------------------------------------------------------------%
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);
numIons  = ORG_STRUC.numIons;

%.LATTICE
para = zeros(3,3);
pos = [ 0, 0, 0];


atomSymbol = cell(1,sumIons);
iC=1;
for typeA = 1 : length(numIons)
    for i = 1 : numIons(typeA)
        atomSymbol{iC}=megaDoof( ORG_STRUC.atomType(typeA) );
        iC = iC + 1;
    end
end



%-----------------------------------------------------
cd(ORG_STRUC.homePath);
fpath = [ORG_STRUC.resFolder '/PATH/path' num2str(POP_STRUC.step) '.xyz'];
fp = fopen(fpath, 'w');


%-- ASXF HEAD Format

for i = 1:numImages
    fprintf(fp,' %3d\n', sumIons);
    fprintf(fp, ' Image-%2d ,\n', i);
    for j = 1:sumIons
        fprintf(fp, ' %2s   %10.6f   %10.6f  %10.6f\n', atomSymbol{j}, POP_STRUC.POPULATION(i).CARTECOORDS(j,:) );
    end
end

fclose(fp);
