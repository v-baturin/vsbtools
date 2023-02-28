function ANS = Random_200(Ind_No, Vacuum, cell)

global ORG_STRUC
global OFF_STRUC

thicknessS     = ORG_STRUC.thicknessS;
atomType       = ORG_STRUC.atomType;
minDistMatrice = ORG_STRUC.minDistMatrice;

%---------Step1: to find the cell size and lattice
bulk_lat    =ORG_STRUC.bulk_lat;
bulk_pos    =ORG_STRUC.bulk_pos;
bulk_numIons=ORG_STRUC.bulk_ntyp;

[bulk_pos, bulk_lat, bulk_numIons, Trans]=...
Make_MultiCell(bulk_pos, bulk_lat, bulk_numIons, cell);

bulk_atyp = [];
for i = 1:length(bulk_numIons)
    bulk_atyp = [bulk_atyp; atomType(i)*ones(bulk_numIons(i),1)];
end

sur_lat_input = bulk_lat;
sur_lat_input(3,3) = thicknessS;
sur_numIons  = ORG_STRUC.numIons*det(Trans);

%%%%Step 2: to get the atomic positions
goodBad = 0;
Groups = ORG_STRUC.nsym;

while ~goodBad 
     tmp = find(Groups > 0);
     nsym = tmp(ceil(rand*length(tmp)));
     
     cd([ORG_STRUC.homePath '/CalcFoldTemp'])
     [sur_candidate, sur_lat, errorS] = ...
     symope_GB(nsym, sur_numIons, sur_lat_input, minDistMatrice);
     cd(ORG_STRUC.homePath)
     
     if errorS == 0
        [lat, candidate, numIons, chanAList] = ...
        makeSurface(sur_lat, sur_candidate, sur_numIons, bulk_lat, bulk_pos, bulk_numIons, Vacuum);
        goodBad = distanceCheck(candidate, lat, numIons, minDistMatrice);
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
