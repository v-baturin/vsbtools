function content = POSCARContent_Surface(atomType, count, symg, numIons, lattice, coor, lat_bulk, thicknessB, step)


% rewrite the POSCAR writor 

symg='no SG';

if size(coor,1) ~= sum(numIons)
   USPEXmessage(530,'',0);
   disp(pwd);
end

lat_6 = latConverter(lattice);  %1*6
lat_6(4:6) = lat_6(4:6)*180/pi;

header = sprintf('EA%-4d %6.3f %6.3f %6.3f %5.2f %5.2f %5.2f Sym.group: %4d', count, lat_6(1:6), symg);

content = {
    header,
    '1.0'
    };

for latticeLoop = 1 : 3
   content(end+1) = { sprintf('%12.6f %12.6f %12.6f', lattice(latticeLoop,:)) };
end

atomTypeStr = [];
if sum(atomType)>0
   for i=1:length(numIons)
       if numIons(i) > 0
          atomTypeStr =  [atomTypeStr, sprintf('%4s', megaDoof(ceil(atomType(i))))];
       end
   end
     content(end+1) = { atomTypeStr };
end

content(end+1) = {'Selective dynamics'};
content(end+1) = {'Direct'};
for coordLoop = 1 : sum(numIons)
    if chanAList(coordLoop) == 1
       content(end+1) = { sprintf('%12.6f %12.6f %12.6f  T  T  T', coor(coordLoop,:)) };
    else
       if nargin == 9
          content(end+1) = { sprintf('%12.6f %12.6f %12.6f  F  F  F', coor(coordLoop,:)) };
       else
          if coor(i,3)*lat_6(3) > lat_bulk(3) - thicknessB
             content(end+1) = { sprintf('%12.6f %12.6f %12.6f  T  T  T', coor(coordLoop,:)) };
          else
             content(end+1) = { sprintf('%12.6f %12.6f %12.6f  F  F  F', coor(coordLoop,:)) };
          end
       end
    end
end




