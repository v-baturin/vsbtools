function [comp, enth, volume] = parseIndividuals(id_list)
% $Rev$
% $Author$
% $Date$

global ORG_STRUC
global USPEX_STRUC

comp   = [];
enth   = [];
volume = [];
for i=1:length(id_list)
    if ORG_STRUC.molecule == 1
        comp   = [comp; USPEX_STRUC.POPULATION(i).numMols];
    else
        comp   = [comp; USPEX_STRUC.POPULATION(i).numIons];
    end
    enth   = [enth  ; USPEX_STRUC.POPULATION(i).Enthalpies(end)];
    volume = [volume; USPEX_STRUC.POPULATION(i).Vol];
end

end
