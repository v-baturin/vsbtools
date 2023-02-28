function [ Coord ] = getAtomCoords( raw_coords )
% This function is to read the atom coordinates from a raw_coords
% ex: raw_coords = '0.1567(2)' Coord = 0.1567
 Coord = 0;
 coords = '';
    for i=1:length(raw_coords)
       if (raw_coords(i) ~='(')
           coords = strcat(coords, raw_coords(i));
       else if raw_coords(i)=='('
               break;
           end
             
       end
    end
    
    Coord= str2num(coords);
end

