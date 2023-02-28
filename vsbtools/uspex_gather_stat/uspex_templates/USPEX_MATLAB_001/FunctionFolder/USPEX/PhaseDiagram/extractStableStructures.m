function extractStableStructures(ids)
% $Rev$
% $Author$
% $Date$

global ORG_STRUC
global USPEX_STRUC

% Save found structures to stablePOSCARS_pressure:
unixCmd('rm -f POSCAR stablePOSCARS_pressure');

for i=1:length(ids)
    Write_POSCAR(ORG_STRUC.atomType, ids(i),             ...
                 USPEX_STRUC.POPULATION(ids(i)).symg,    ...
                 USPEX_STRUC.POPULATION(ids(i)).numIons, ...
                 USPEX_STRUC.POPULATION(ids(i)).LATTICE, ...
                 USPEX_STRUC.POPULATION(ids(i)).COORDINATES);

    unixCmd('cat POSCAR >> stablePOSCARS_pressure');
end

unixCmd('rm -f POSCAR');

end
