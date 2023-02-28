function POP = Random_Init_311(Ind_No, numMols)


global ORG_STRUC

CminD   = ORG_STRUC.CenterminDistMatrice;
 minD   = ORG_STRUC.minDistMatrice;

%-----------------------------LATTICE------------------------------
lat1 = numMols*ORG_STRUC.latVolume';
[typesAList, MtypeLIST, numIons] = GetPOP_MOL(numMols);
%-----------------------------LATTICE------------------------------
 
%--------------------------Initilization---------------------------
goodStruc     = 0;
badSymmetry   = 0;
Freq_max      = 5;
tmp = find(ORG_STRUC.nsym > 0);
nsym = tmp(ceil(rand*length(tmp))); % pick a random group 
%--------------------------Initilization----------------------------
tic
time  = 100;  %If time exceeds 100 second, need to change parameters
time0 = 100;

%------------------------Generation Loop ---------------------------
while ~goodStruc

      if badSymmetry > Freq_max
          tmp = find(ORG_STRUC.nsym > 0);
          nsym = tmp(ceil(rand*length(tmp))); % pick a random group 
          badSymmetry = 0;
          if (toc > time0)
             fac = sqrt(toc/time);
             time = toc;
             disp(['']);
             disp(['Long time to generate structure: '...
                    num2str(toc, '%4d') ' seconds']);
             CminD = CminD*fac;
             fac3  = fac*fac*fac;
             lat1  = lat1*fac3;
             disp(['Increasing Molcenter by ' num2str(fac) ' times']);
             disp(['=====>      ' num2str(CminD, '%8.3f')]);
             disp(['Increasing volume by ' num2str(fac3) ' times']);
             disp(['=====>      ' num2str(lat1,   '%8.3f')]);
             ORG_STRUC.latVolume = ORG_STRUC.latVolume*fac3;
             disp(['']);
             ORG_STRUC.CenterminDistMatrice = CminD;
          end

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
                    POP.MOLECULES  = Molecules;
                    POP.numMols    = numMols;
                    POP.MtypeLIST  = MtypeLIST;
                    POP.typesAList = typesAList;
                    POP.numIons    = numIons;
                    POP.LATTICE    = lat;
                    POP.howCome    = '  Random  ';

                    disp( ['Structure '  num2str(Ind_No) ...
                    ' built (Z = ' num2str(numMols) ') with the symmetry '... 
                                num2str(nsym) ' (' spaceGroups(nsym)  ') ']);
                      
                    goodStruc = 1;
                    break
                    
                end
             end
          end
      end
end
