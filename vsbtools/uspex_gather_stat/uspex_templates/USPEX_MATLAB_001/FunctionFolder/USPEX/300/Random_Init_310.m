function [ newCoords, lat] = Random_Init_310(Ind_No, numMols)

% implemented - USPEX Version 9.3.7
global ORG_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CminD   = ORG_STRUC.CenterminDistMatrice;
 minD   = ORG_STRUC.minDistMatrice;
atomType= ORG_STRUC.atomType;
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
Freq_max      = 5;
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

   
      cd(['CalcFoldTemp']);
      [candidate, lat, numSites, Operation, errorMsg] = ...
                     symope_3D_MOL(nsym, numMols, lat1, CminD);
      cd(ORG_STRUC.homePath);
      if (errorMsg == 2) && (ORG_STRUC.minAt < ORG_STRUC.maxAt)
          ORG_STRUC.nsym(nsym)=0; 
          badSymmetry = Freq_max + 1;
          
      elseif errorMsg == 0
          if distanceCheck(candidate, lat, numMols, CminD-0.2)
             for item=1:20

                Molecules = GetOrientation(candidate, lat, numSites,...
                                           Operation, MtypeLIST, nsym);

                if newMolCheck(Molecules,lat, MtypeLIST, minD);
                    disp( ['Structure '  num2str(Ind_No) ...
                    ' built (Z = ' num2str(numMols) ') with the symmetry '... 
                                num2str(nsym) ' (' spaceGroups(nsym)  ') ']);
                      
                    newCoords = GetCoor(typesAList, Molecules, lat, atomType);
                    goodStruc = 1;
                    break
                    
                end
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

