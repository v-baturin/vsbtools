function [candidate_2D, lat_2D, candidate, lat] = Random_Init_M200(Ind_No, numIons)

% implemented - USPEX Version 8.5.0
global ORG_STRUC

newSym = 1;
badSymmetry = 0;
goodBad = 0;
minD        = ORG_STRUC.minDistMatrice;
thickness   = ORG_STRUC.thicknessS;

if ORG_STRUC.constLattice
    if ORG_STRUC.minAt ~= ORG_STRUC.maxAt
        disp('Single block does not support constantLattice');
        disp('Please check your INPUT.txt file, quit......');
        quit;
    end
    
else
    area = ORG_STRUC.latVolume/(thickness+2.0);
    area = area*sum(numIons)/sum(ORG_STRUC.numIons);
end

while ~goodBad
    if (badSymmetry > 15)
        badSymmetry = 0;
        newSym = 1;      % change the symmetry group if can't generate the crystal
    end
    
    badSymmetry = badSymmetry + 1;
    
    if newSym
        tmp = find(ORG_STRUC.nsym > 0);
        nsym = tmp(ceil(rand*length(tmp))); % pick a random group
        newSym = 0;
    end
    
    cd([ORG_STRUC.homePath '/CalcFoldTemp'])
    if ORG_STRUC.constLattice
        lat  = ORG_STRUC.lattice;
        [candidate_2D, lat_2D, errorS] = symope_GB(nsym, numIons, lat, minD);
    else
        [candidate_2D, lat_2D, errorS] = symope_2D(nsym, numIons, area, minD, thickness);
    end
    cd(ORG_STRUC.homePath)
    
    if errorS == 0
        [lat,candidate] = make2D(lat_2D, candidate_2D, ORG_STRUC.vacuumSize(1));
        goodBad = distanceCheck(candidate, lat, numIons, minD);
        if goodBad == 1
            goodBad = checkConnectivity(candidate, lat, numIons);
        end
    end
    
    if goodBad
        disp(['Structure ' num2str(Ind_No) ' generated with plane group: '...
            num2str(nsym)]);
        if ORG_STRUC.minAt < ORG_STRUC.maxAt
            disp(['numIons: ' num2str(numIons)]);
        end
    end
end