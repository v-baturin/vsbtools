function [candidate_2D, lat_2D, candidate, lat] = Random_Init_M210(Ind_No, numMols)

global ORG_STRUC

CminD      = ORG_STRUC.CenterminDistMatrice-1;
minD       = ORG_STRUC.minDistMatrice;
vacuum     = ORG_STRUC.vacuumSize(1);
atomType   = ORG_STRUC.atomType;
thickness2 = ORG_STRUC.thicknessS;
thickness1 = ORG_STRUC.thicknessB;

%-----------------------------LATTICE------------------------------
if ORG_STRUC.constLattice    % fixed lattice
    lat1 = ORG_STRUC.lattice;
else                         % volume
    lat1 = numMols*ORG_STRUC.latVolume;
end
[typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
%-----------------------------LATTICE------------------------------

%--------------------------Initilization---------------------------
goodStruc     = 0;
badSymmetry   = 0;
Freq_max      = 15;
tmp = find(ORG_STRUC.nsym > 0);
nsym = tmp(ceil(rand*length(tmp))); % pick a random group 
%--------------------------Initilization----------------------------

while ~goodStruc

    if badSymmetry > Freq_max
       tmp = find(ORG_STRUC.nsym > 0);
       nsym = tmp(ceil(rand*length(tmp))); % pick a random group 
       badSymmetry = 0;
    else
       badSymmetry = badSymmetry + 1;
    end

    
     cd([ORG_STRUC.homePath '/CalcFoldTemp'])
     [coor, lat_2D, numSites, Operation, errorMsg] = ...
     symope_2D_MOL(nsym, numMols, lat1, CminD+1, thickness1);
     cd(ORG_STRUC.homePath)

     if errorMsg == 0
        [lat_2D,coor] = make2D(lat_2D, coor, thickness2);              
        for item=1:30
            Molecules = GetOrientation(coor, lat_2D, numSites, ...
                                        Operation, MtypeLIST, nsym);
            if ~newMolCheck(Molecules,lat_2D, MtypeLIST, minD+0.4);  
                disp(['Structure ' num2str(Ind_No) ...
                 ' generated from Random symmetry ' num2str(nsym) ]);
                
                candidate_2D = GetCoor(typesAList, Molecules, lat_2D, atomType);
                [lat,candidate] = make2D(lat_2D, candidate_2D, vacuum);         
                goodStruc = 1;
                break
            end
        end
     end 
end 
%---------------------------------------------------------------
function newCoords = GetCoor(typesAList, Molecules, lat, atomType)

cattedCoors = [];
for j = 1: length(Molecules)
    cattedCoors = cat(1,cattedCoors,Molecules(j).MOLCOORS);
end
saveded = cattedCoors/lat;
newCoords = zeros(0,3);
for m = 1: length(atomType)
    s = find(typesAList == atomType(m));
    newCoords = cat(1,newCoords, saveded(s,:));
end

