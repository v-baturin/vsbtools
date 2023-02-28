function [lat,candidate,errorS] = random_symmetry(Ind_No,latVolume,numIons,nsym,minDistMatrice)

global ORG_STRUC
global POP_STRUC

if ORG_STRUC.constLattice    % lat1 = lattice
    if ORG_STRUC.minAt ~= ORG_STRUC.maxAt
        error('Single block or varcomp doesnot support constantLattice, Please check your INPUT.txt file')
        quit();
    end
    lat1 = ORG_STRUC.lattice;
    lat = lat1;
else                         % lat1 = lattice volume
    lat1 = latVolume;
    new = 0;
    while ~new
        startLat = rand(6,1);
        startLat(4:6) = startLat(4:6)*(pi/2);
        check_startLat = latConverter(startLat);
        volLat = det(check_startLat);
        if volLat == 0
            continue;
        end
        % scale the lattice to the volume we asume it approximately to be
        ratio = lat1/volLat;
        startLat(1:3) = startLat(1:3)*(ratio)^(1/3);
        [dummy, lat2] = optLattice([0 0 0], latConverter(startLat)); % optimize lattice
        new = latticeCheck(lat2);
        if new
            lat = lat2;
        end
    end
end
% H. Stokes code to create a crystal with given symmetry
% OLD: if nsym == 1 - just keep the random candidate that we generated
if sum(ORG_STRUC.splitInto) > 3 %splitcell
    [lat, errorS, candidate] = splitBigCell(latConverter(lat), ...
    ORG_STRUC.splitInto(ceil(length(ORG_STRUC.splitInto)*rand)),numIons, nsym);
else
    cd([ORG_STRUC.homePath '/CalcFoldTemp']);
    [candidate, lat, errorS] = symope_crystal(nsym, numIons, lat1, minDistMatrice, ORG_STRUC.sym_coef);
    cd(ORG_STRUC.homePath)
end
%if errorS == 0
%    cd([ORG_STRUC.homePath '/CalcFoldTemp']);
%    [resulted_sg_no, nothing] = determine_spacegroup(lat, candidate, ORG_STRUC.atomType, ...
%        numIons, nsym, POP_STRUC.bodyCount + Ind_No, ...
%        ORG_STRUC.SGtolerance, ORG_STRUC.homePath);
%    cd(ORG_STRUC.homePath)
%else
%    resulted_sg_no = 0;
%end

