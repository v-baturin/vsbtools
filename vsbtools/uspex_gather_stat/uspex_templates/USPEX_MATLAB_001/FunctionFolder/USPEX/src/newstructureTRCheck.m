function [newstructure] = newstructureTRCheck(lat,candidate)

%store generated structures and drop secondly apeared ones
%USPEX_STRUC.EXT_RandTop.STRUCTURES is global structure containing information
%on structures generated with topological random

global USPEX_STRUC

newstructure = true;
for i = 1:length(USPEX_STRUC.EXT_RandTop.STRUCTURES)
    l1 = USPEX_STRUC.EXT_RandTop.STRUCTURES(i).LATTICE;
    l2 = lat;
    if (norm(l1-l2) < 0.05*max(norm(l1),norm(l2))+1)&&isequal(USPEX_STRUC.EXT_RandTop.STRUCTURES(i).COORDINATES,candidate)
        newstructure = false;
        break;
    end
end
if newstructure
    USPEX_STRUC.EXT_RandTop.STRUCTURES(end+1).LATTICE = lat;
    USPEX_STRUC.EXT_RandTop.STRUCTURES(end).COORDINATES = candidate;
end



%if length(USPEX_STRUC.EXT_RandTop.STRUCTURES) == 0
%    USPEX_STRUC.EXT_RandTop.STRUCTURES(end+1).LATTICE = lat;
%    USPEX_STRUC.EXT_RandTop.STRUCTURES(end).COORDINATES = candidate;
%else
%    for i = 1:length(USPEX_STRUC.EXT_RandTop.STRUCTURES)
%        l1 = USPEX_STRUC.EXT_RandTop.STRUCTURES(i).LATTICE;
%        l2 = lat;
%        if (norm(l1-l2) < 0.05*max(norm(l1),norm(l2))+1)&&isequal(USPEX_STRUC.EXT_RandTop.STRUCTURES(i).COORDINATES,candidate)
%            newstructure = false;
%            break;
%        end
%    end
%    if newstructure
%        USPEX_STRUC.EXT_RandTop.STRUCTURES(end+1).LATTICE = lat;
%        USPEX_STRUC.EXT_RandTop.STRUCTURES(end).COORDINATES = candidate;
%    end
%end

