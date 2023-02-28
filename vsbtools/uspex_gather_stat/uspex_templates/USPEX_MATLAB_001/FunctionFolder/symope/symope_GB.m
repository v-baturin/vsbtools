function [candidate, lattice, errorS] = symope_GB(nsym, numIons, lattice, minD)

global ORG_STRUC

spgBINDIR  =[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];
fixRndSeed = ORG_STRUC.fixRndSeed;
Write_Stokes_input_GB(numIons, minD, nsym, lattice, fixRndSeed);
unixCmd([spgBINDIR  '/random_2d < rc.in > rc.out']);
[candidate, lattice, errorS] = Read_Stokes_output('rc.out');
lattice(4:6) = lattice(4:6)*pi/180;
lattice = latConverter(lattice);


function Write_Stokes_input_GB(numIons, minD, nsym, lattice, fixRndSeed);

To_Delete = find(numIons==0);
numIons(To_Delete) = [];
minD(To_Delete,:) = [];
minD(:,To_Delete) = [];

Lat = latConverter(lattice);
thickness = norm(Lat(3));
area = det(lattice)/thickness;

fp = fopen('rc.in', 'w');
fprintf(fp,'%4d   ! space group \n', nsym);
fprintf(fp,'%6.3f ! area of primitive unit cell\n', area);
fprintf(fp,'%6.3f %6.3f %6.3f ! vector of primitive unit cell\n', ...
                                       Lat(1), Lat(2), Lat(6)*180/pi);
fprintf(fp,'%6.3f ! thickness of layer\n', thickness-1);
fprintf(fp,'%4d   ! number of types of atoms\n', length(numIons));

for i=1:length(numIons)
    fprintf(fp,'%4d ', numIons(i));
end

fprintf(fp, '! number of atoms of each type\n');

for ii = 1 : size(minD,1)
    for jj = 1 : size(minD,2)
        fprintf(fp, '%5.3f ', minD(ii,jj));
    end
end
fprintf(fp, '! minimum distance between atoms \n');
fprintf(fp, '0  ! symmetry operation \n');
if fixRndSeed>0
    fprintf(fp, '%10d %10d ! RandSeeds \n', round(rand(1,2)*10^6));
end
fclose(fp);

