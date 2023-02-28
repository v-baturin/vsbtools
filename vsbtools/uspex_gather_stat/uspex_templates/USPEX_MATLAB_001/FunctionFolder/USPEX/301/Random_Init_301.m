function [candidate, lat] = Random_Init_301(Ind_No, numBlocks, topologyRandom)

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CREATING random structures using space groups provided  or toplogy %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodStructure = 0;
goodLattice = 0;

newSym = 1;
badSymmetry = 0;
failedDist = 0;
minDistMatrice = ORG_STRUC.minDistMatrice;
badStruc = 0;
badLat   = 0;

numIons   = numBlocks*ORG_STRUC.numIons;

tic
while goodStructure + goodLattice  ~= 2

    failedTime = toc;
    if (failedDist > 10000) || (failedTime > 300)
        if minDistMatrice(1,1) > 0.8*ORG_STRUC.minDistMatrice(1,1)
            if failedTime > 300
                USPEXmessage(501,'',0);
            else
                USPEXmessage(502,'',0);
            end
            %disp('Please check your IonDistances parameter.');
            %disp(' ');
            minDistMatrice = 0.9*minDistMatrice;
            failedDist = 0;
            tic
        else
            USPEXmessage(562,'',0);
            quit;
        end
    end

    errorS = 0;
    
    latVolume = dot(numBlocks,ORG_STRUC.latVolume);

    if topologyRandom == 0    
        if (badSymmetry > 15) || ((badSymmetry > 5) && (sum(ORG_STRUC.splitInto)>3))
            badSymmetry = 0;
            newSym = 1;      % change the symmetry group if can't generate the crystal
        end
        badSymmetry = badSymmetry + 1;
        if newSym
            tmp = find(ORG_STRUC.nsym > 0);
            nsym = tmp(ceil(rand*length(tmp))); % pick a random group from those specified by user
            newSym = 0;
        end
        [lat,candidate,errorS] = random_symmetry(Ind_No,latVolume,numIons,nsym,minDistMatrice);
    else
        %get coordination numbers from interpolation grid
        %coordinationNumbers = coordinationNumbersForComposition(numIons);
        %disp(str_coordinationNumbers);
        %disp(str_numIons);
        %disp(str_latVolume);
        newstructure = false;
        top_count = 0;
        elapsedTime = toc;
        while (~newstructure)&&(top_count < 100)&&(elapsedTime<120)
            [topology_number, candidate, lat] = random_topology(latVolume,numIons);
            newstructure = distanceCheck(candidate, lat, numIons, minDistMatrice)&latticeCheck(lat)&newstructureTRCheck(lat,candidate);
            top_count = top_count + 1;
            elapsedTime = toc;
        end
        if ~newstructure
            topologyRandom = 0;
            tic;
            errorS = 1;
            %lat = zeros(3,3);
            %candidate = zeros(sum(numIons),3);
        end
    end


    if errorS == 0
        goodStructure = distanceCheck(candidate, lat, numIons, minDistMatrice*ORG_STRUC.sym_coef);
        goodLattice = latticeCheck(lat);
        badLat = badLat + 1 - goodLattice;
        badStruc = badStruc + 1 - goodStructure;
        if badLat > 1000
          USPEXmessage(503,'',0);
          badLat = 0;
        end
        if badStruc > 1000
            USPEXmessage(504,'',0);
            badStruc = 0;
        end
    else
        goodStructure = 0;
        goodLattice = 0;
    end

    if goodStructure + goodLattice == 2
        if topologyRandom == 0
            disp(['Structure ' num2str(Ind_No) ' built with the symmetry group ' num2str(nsym) ' (' spaceGroups(nsym) ') , composition ' num2str(numIons)]);
        else
            disp(['Structure ' num2str(Ind_No) ' built with the topology ' num2str(topology_number) ' , composition ' num2str(numIons) ]);
        end
    else
        failedDist = failedDist + 1;
    end

end

