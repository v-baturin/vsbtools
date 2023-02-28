function Write_Stokes_input(numIons, minD, nsym, lat1, fixRndSeed, sym_coef);

To_Delete = find(numIons==0);
numIons(To_Delete) = [];
minD(To_Delete,:) = [];
minD(:,To_Delete) = [];

fp = fopen('rc.in', 'w');
fprintf(fp,'%4d ! space group \n', nsym);
fprintf(fp,'%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f ! lattice of primitive unit cell\n', lat1(1:6));
fprintf(fp,'%4d ! number of types of atoms\n', length(numIons));

for i=1:length(numIons)
    fprintf(fp,'%4d ', numIons(i));
end
fprintf(fp, '! number of atoms of each type\n');

if ~isempty(sym_coef)
   minD = minD*sym_coef;
end
for ii = 1 : size(minD,1)
    for jj = 1 : size(minD,2)
        fprintf(fp, '%5.3f ', minD(ii,jj));
    end
end
fprintf(fp, '! minimum distance between atoms \n');

if ~isempty(sym_coef)
   fprintf(fp, ' 1 coefficient between minDist and symmetrization distance \n');
end

if fixRndSeed>0
    fprintf(fp, '%10d %10d ! RandSeeds \n', round(rand(1,2)*10^6));
end
fclose(fp);
