function [newlat, newcoord, numIons, operation, Error] = find_symmetry(lat, coord, numIons, atomType, tolerance)

global ORG_STRUC
USPEXPath = ORG_STRUC.homePath;
getPy    =[ORG_STRUC.USPEXPath '/FunctionFolder/getInput.py'];
spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup'];
newlat = lat;
newcoord = coord;
operation = [];
Error = 0;

%Step1: prepare the input
Write_Symmetry_input(lat, coord, numIons, atomType, tolerance);
[a, b] = unix([spgBINDIR '/findsym_new < sym.in > sym.out']);

%Step2: obtain the lattice and wycoff position, symmetry operations;
[lat_6, REF_COOR, operation, Error] = Read_CIF('sym.out');
if ~Error
   newlat = latConverter(lat_6');
   N_wyckoff = size(REF_COOR,1);
   Num = size(operation,1)/4;    %Number of operations
   for i = 1:Num
       Ope = operation(i*4-3:i*4,:);
       for k = 0:N_wyckoff-1
           newcoord(i+k*Num,:) = REF_COOR(k+1,:)*(Ope(1:3,:)) + Ope(4,:);
       end
   end
   numIons = Num*N_wyckoff;
end

newcoord = RemoveAtom(newcoord);

function coor1 = RemoveAtom(coor)

coor = coor - floor(coor);
coor1 = coor(1,:);

for i = 2:size(coor,1)
    Add = 1;
    for j = 1:size(coor1,1)
        if norm(coor(i,:) - coor1(j,:)) < 0.05
           Add = 0;
           break;
        end
    end
    if Add == 1
       coor1 = [coor1; coor(i,:)];
    end
end
