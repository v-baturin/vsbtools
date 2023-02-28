function [candidate, lattice, numSites, Operation, errorS] = symope_2D_MOL(nsym, numIons, volume, minD, thickness)

global ORG_STRUC

spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];
fixRndSeed = ORG_STRUC.fixRndSeed;

Write_Stokes_input_2D_MOL(numIons, minD, nsym, volume, thickness, fixRndSeed);
unixCmd([spgBINDIR  '/random_2d < rc.in > rc.out']);
[Operation, numSites, candidate, lattice, errorS] = Read_Stokes_output_MOL('rc.out');
if ~errorS
   lattice(4:6) = lattice(4:6)*pi/180;
   lattice = latConverter(lattice);
end


function Write_Stokes_input_2D_MOL(numIons, minD, nsym, volume, thickness, fixRndSeed)

To_Delete = find(numIons==0);
numIons(To_Delete) = [];
minD(To_Delete,:) = [];
minD(:,To_Delete) = [];

if size(volume) == [3,3];
  Lat_tmp = latConverter(volume);
  lat_in = [Lat_tmp(1), Lat_tmp(2), Lat_tmp(4)*180/pi];
  volume = det(volume)/thickness;
else
  lat_in = [1 1 90];
  volume = volume/thickness;
end

fp = fopen('rc.in', 'w');
fprintf(fp,'%4d   ! space group \n', nsym);
fprintf(fp,'%6.3f ! area of primitive unit cell\n', volume);
fprintf(fp,'%6.3f %6.3f %6.3f ! vector of primitive unit cell\n', lat_in);
fprintf(fp,'%6.3f ! thickness of layer\n', thickness);
fprintf(fp,'%4d   ! number of types of atoms\n', length(numIons));

for i=1:length(numIons)
   fprintf(fp,'%4d ', numIons(i));
end

fprintf(fp, '! number of atoms of each type\n');

minD = minD;
for ii = 1 : size(minD,1)
 for jj = 1 : size(minD,2)
  fprintf(fp, '%5.3f ', minD(ii,jj));
 end
end
fprintf(fp, '! minimum distance between atoms \n');
fprintf(fp, '1    ! symmetry operation \n');
if fixRndSeed>0
    fprintf(fp, '%10d %10d ! RandSeeds \n', round(rand(1,2)*10^6));
end
fclose(fp);

