function [candidate, lat] = Random_Init_300(Ind_No, numIons, topologyRandom)

% implemented - USPEX Version 8.5.0

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% CREATING random structures using space groups provided or topology %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
goodStructure = 0;
goodLattice = 0;

newSym = 1;
badSymmetry = 0;
failedDist = 0;
minDistMatrice = ORG_STRUC.minDistMatrice;
badStruc = 0;
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
            minDistMatrice = 0.9*minDistMatrice;
            failedDist = 0;
            tic
        else
            USPEXmessage(562,'',0);
            quit;
        end
    end
    
    errorS = 0;

    latVolume = ORG_STRUC.latVolume*sum(numIons)/sum(ORG_STRUC.numIons);

    if topologyRandom == 0
        if (badSymmetry > 15) || ((badSymmetry > 5) && (sum(ORG_STRUC.splitInto)>3))
            badSymmetry = 0;
            newSym = 1;      % change space group if can't generate the crystal
        end
        badSymmetry = badSymmetry + 1;
        if newSym
            tmp = find(ORG_STRUC.nsym > 0);
            nsym = tmp(ceil(rand*length(tmp))); % pick a random group from INPUT
            newSym = 0;
        end
        [lat, candidate, errorS] = random_symmetry(Ind_No,latVolume,numIons,nsym,minDistMatrice);
    else
        %get coordination numbers from interpolation grid
        %coordinationNumbers = coordinationNumbersForComposition(numIons);
        %disp(str_coordinationNumbers);
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
        
        badStruc = badStruc + 1 - goodStructure;
        
        if badStruc > 1000
            %disp('Failed 1000 times during random structure generation.')
            %disp(' ')
            USPEXmessage(503,'',0);
            badStruc = 0;
        end
    else
        goodStructure = 0;
        goodLattice = 0;
    end
    
    if goodStructure + goodLattice == 2
        %if topologyRandom == 0
        %    % Structures with wrong space group:
        %    %if nsym >= 1
        %    %   if nsym ~= resulted_sg_no
        %    %      if ~isfield(ORG_STRUC, 'wrong_spacegroups')
        %    %        ORG_STRUC.wrong_spacegroups = 1;
        %    %      else
        %    %        ORG_STRUC.wrong_spacegroups = ORG_STRUC.wrong_spacegroups + 1;
        %    %      end
        %    %   end
        %    %end
        %end

        if ORG_STRUC.maxAt > ORG_STRUC.minAt
            if topologyRandom == 0
                disp(['Structure ' num2str(Ind_No) ' built with space group '...
                num2str(nsym) ' (' spaceGroups(nsym) ...
                ') , composition ' num2str(numIons)]);
            else
                disp(['Structure ' num2str(Ind_No) ' built with topology '...
                num2str(topology_number) ...
                ' , composition ' num2str(numIons) ]);
            end
        else
            if topologyRandom == 0
                disp(['Structure ' num2str(Ind_No) ' built with space group '...
                num2str(nsym) ' (' spaceGroups(nsym) ')']);
            else
                disp(['Structure ' num2str(Ind_No) ' built with topology '...
                num2str(topology_number)  ]);
            end
        end
    else
        failedDist = failedDist + 1;
    end
    
end
