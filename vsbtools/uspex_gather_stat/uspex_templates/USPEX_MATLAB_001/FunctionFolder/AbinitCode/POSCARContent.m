function content = POSCARContent(atomType, count, symg, numIons, lattice, coor, order)


% This function is designed for store the POSCAR format cell, easy for POCAR file write
%  
% count : 1. when count is a string, the string will be used as the 1st line of POSCAR 
%         2. when count is a number, 'EA%-4d %6.3f %6.3f %6.3f %5.2f %5.2f %5.2f' format 
%            will be used as the 1st line
% atomType:  0: do not write element symbol (VASP-4.6 format)
%        array:        write element symbol (VASP-5 format)

if size(coor,1) ~= sum(numIons)
   USPEXmessage(530,'',0);
   disp(pwd);
end

if isempty(symg)
   symg = 0;
end

lat_6 = latConverter(lattice);  %1*6
lat_6(4:6) = lat_6(4:6)*180/pi;

if ischar(count)
  header = count;
else
  header = sprintf('EA%-4d %6.3f %6.3f %6.3f %5.2f %5.2f %5.2f Sym.group: %4d', count, lat_6(1:6), symg);
end


content = {
    header,
    '1.0'
    };

for latticeLoop = 1 : 3
   content(end+1) = { sprintf('%12.6f %12.6f %12.6f', lattice(latticeLoop,:)) };
end

atomTypeStr = '';
if sum(atomType) > 0  %QZ: only write when atomType ~=0
   %atomType(find(atomType==0))=[];
   for i=1:length(atomType)
        if numIons(i)>0
          atomTypeStr =  [atomTypeStr, sprintf('%4s', megaDoof(ceil(atomType(i))))];
        end
   end
   content(end+1) = { atomTypeStr };
end

numIonsStr  = '';
%if sum(numIons) > 0
   %numIons(find(numIons==0))=[];
   for i=1:length(numIons)
        if numIons(i)>0
          numIonsStr  =  [numIonsStr,  sprintf('%4d', numIons(i)) ];
        end
   end
   content(end+1) = { numIonsStr };
%end

content(end+1) = {'Direct'};
for coordLoop = 1 : sum(numIons)
    if nargin == 7
       content(end+1) = { sprintf('%12.6f %12.6f %12.6f', coor(coordLoop,:), order(coordLoop)) };
    else
       content(end+1) = { sprintf('%12.6f %12.6f %12.6f %12.6f', coor(coordLoop,:)) };
    end
end




