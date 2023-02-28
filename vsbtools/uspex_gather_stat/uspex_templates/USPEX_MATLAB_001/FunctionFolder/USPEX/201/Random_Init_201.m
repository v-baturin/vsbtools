function ANS = Random_201(Ind_No, Vacuum, cell)

global ORG_STRUC
global OFF_STRUC
thicknessS     = ORG_STRUC.thicknessS;
atomType       = ORG_STRUC.atomType;
minDistMatrice = ORG_STRUC.minDistMatrice;
%%%%Step1: to find the cell size and lattice
bulk_lat    =ORG_STRUC.bulk_lat;
bulk_pos    =ORG_STRUC.bulk_pos;
bulk_atyp   =ORG_STRUC.bulk_atyp;
bulk_numIons=ORG_STRUC.bulk_ntyp;

[bulk_pos, bulk_lat, bulk_numIons, Trans]=...
Make_MultiCell(bulk_pos, bulk_lat, bulk_numIons, cell);

count = 0;
for i = 1:length(bulk_numIons)
   bulk_atyp(count+1:count+bulk_numIons(i)) = atomType(i);
   count = count + bulk_numIons(i);
end

sur_lat_input = bulk_lat;
sur_lat_input(3,3) = thicknessS;

sur_numIons=[];
N_T = length(atomType);
for i=1:N_T
    if (ORG_STRUC.numIons(i)==0)
        Surface_numIons(i)=0;
    else
        atom = ORG_STRUC.numIons(i)*det(Trans);
        sur_numIons(i)=RandInt(1,1,[1,atom]);
    end
end

%%%%Step 2: to get the atomic positions
goodBad = 0;
Groups = ORG_STRUC.nsym;

while ~goodBad
    tmp = find(Groups > 0);
    nsym = tmp(ceil(rand*length(tmp)));
    cd([ORG_STRUC.homePath '/CalcFoldTemp']);
    [sur_candidate, sur_lat, errorS] = ...
    symope_GB(nsym, sur_numIons, sur_lat_input, minDistMatrice);
    cd(ORG_STRUC.homePath)
    if errorS == 0
       [lat, candidate, numIons, chanAList] = ...
       makeSurface(sur_lat, sur_candidate, sur_numIons, bulk_lat, bulk_pos, bulk_numIons, Vacuum);
       goodBad = distanceCheck(candidate, lat, numIons, ORG_STRUC.minDistMatrice);
       %if ~goodBad
       %   bulk_lat
       %   bulk_numIons
       %   bulk_pos
       %   quit
       %end
    else
       Groups(nsym)=0;
    end
    if goodBad
       disp(['Structure ' num2str(Ind_No) ' generated randomly']);
       disp(['composition: ' num2str(sur_numIons) ';   cell index: ' num2str(cell)]);
       ANS.candidate     = candidate;
       ANS.numIons       = numIons;
       ANS.lat           = lat;
       ANS.chanAList     = chanAList;
       ANS.sur_lat       = sur_lat;
       ANS.sur_candidate = sur_candidate;
       ANS.sur_numIons   = sur_numIons;
       ANS.bulk_lat      = bulk_lat;
       ANS.bulk_pos      = bulk_pos;
       ANS.bulk_numIons  = bulk_numIons;
       break
    end
end
