function [candidate, lat] = Random_Init_000(Ind_No, numIons)

% implemented - USPEX Version 8.5.0

global ORG_STRUC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sym_coef = ORG_STRUC.sym_coef;
goodStructure = 0;

newSym = 1;
badSymmetry = 0;
failedDist = 0;
minD = ORG_STRUC.minDistMatrice;
vacuumSize = ORG_STRUC.vacuumSize(1);
badStruc = 0;
tic

while goodStructure ~= 1

  failedTime = toc;
  if (failedDist > 10000) || (failedTime > 300)
   if minD(1,1) > 0.8*ORG_STRUC.minDistMatrice(1,1)
     if failedTime > 300
         USPEXmessage(505,'',0);
     else
         USPEXmessage(506,'',0);
     end
     minD = 0.9*minD;
     failedDist = 0;
     tic
   else
      USPEXmessage(562,'',0);
      quit;
   end
  end

  errorS = 0;
  startLat = rand(6,1);
  startLat(4:6) = (pi/2);
  check_startLat = latConverter(startLat);
  volLat = det(check_startLat);
  ratio = ORG_STRUC.latVolume/volLat;
  startLat(1:3) = startLat(1:3)*(ratio)^(1/3); 
  lat = latConverter(startLat);
    

  if badSymmetry > 500    % 100k clusters were generated in ~110 sec
    if sqrt(lat(1,1)*lat(2,2)) < lat(3,3)
      coef = (0.65 + 0.4*rand);
    else
      coef = (0.95 + 0.4*rand);
    end
    lat(3,3) = lat(3,3)*coef;
    lat(1,1) = lat(1,1)/sqrt(coef);
    lat(2,2) = lat(2,2)/sqrt(coef);
  elseif badSymmetry > 150
     badSymmetry = 0;
     newSym = 1;      % change the symmetry if can't generate the cluster
  end
  badSymmetry = badSymmetry + 1;
  if newSym
     tmp = ceil(rand*size(ORG_STRUC.nsymN,1));
     nsym = ORG_STRUC.nsym(ORG_STRUC.nsymN(tmp,1):ORG_STRUC.nsymN(tmp,2));
     newSym = 0;
  end
  
  
  
  
  
  [candidate, lat, errorS] = symope_000(nsym, numIons, lat, minD*sym_coef);

  
  
  
  
  
  if errorS == 0
    goodStructure = distanceCheck(candidate, lat, numIons, minD*sym_coef);
    if goodStructure == 1  %we also need to check the connectivity
       goodStructure = checkConnectivity(candidate, lat, numIons); 
    end
    badStruc = badStruc + 1 - goodStructure;
    if badStruc > 1000
      USPEXmessage(504,'',0);
      disp(' ')
      badStruc = 0;
    end
  else
    goodStructure = 0;
  end

  if goodStructure == 1
      [lat, candidate] = reduce_Cluster(lat, candidate);
      [lat, candidate] = makeCluster(lat, candidate, vacuumSize);
      disp(['Structure ' num2str(Ind_No) ' generated with random symmetry: '...
      num2str(nsym)]);
  else
    failedDist = failedDist + 1;
  end

end
