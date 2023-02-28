function Initialize_POP_STRUC_PSO()

% USPEX Version 9.4.0
% modulization

global ORG_STRUC
global POP_STRUC


POP_STRUC = struct('POPULATION',{}, 'SOFTMODEParents',{}, 'SOFTMUTATED',{}, 'resFolder', {},'generation',{}, ...
    'DoneOrder',{}, 'bodyCount', {}, 'ranking',{},'bad_rank',{}, 'fitness', {}, 'PSO', {}, 'bestPSOstruc', {});

POP_STRUC(1).POPULATION = struct('COORDINATES', {}, 'INIT_COORD', {}, 'LATTICE', {}, 'INIT_LAT', {}, 'INIT_numIons', {}, ...
    'order',{}, 'FINGERPRINT', {}, 'K_POINTS', {},'Step', {}, 'Enthalpies', {}, 'Error',{},'Done',{},'ToDo',{},'Parents',{}, ...
    'struc_entr',{}, 'S_order',{}, 'howCome',{},'JobID',{},'Folder',{}, 'numIons', {}, 'softmode',{}, 'Number',{}, 'symg', {});

POP_STRUC.POPULATION(1) = QuickStart(POP_STRUC.POPULATION);

POP_STRUC(1).PSO = struct('lattice', {}, 'coordinates', {}, 'order', {}, 'numIons', {}, 'fitness', {}, 'fingerprint', {}, 'enthalpy', {}, 'Number', {});
POP_STRUC(1).SOFTMUTATED = struct('FINGERPRINT',{}, 'mutatedAt', {}, 'fallBack', {});
POP_STRUC(1).SOFTMODEParents=struct('lattice',{},'coordinates',{},'fingerprint',{},'eignFre',{},'eignVec',{},'Softmode_Fre',{},'Softmode_num',{},'numIons',{});

POP_STRUC.generation = 1;
POP_STRUC.bodyCount = 0;
POP_STRUC.bad_rank = 0;
POP_STRUC.bestPSOstruc = 1;

%create good initial population. Every individual fulfills hard constraints.
goodPop = 1;
newSym = 1;
badSymmetry = 0;
nsym = 0;
sym_coef = 1;
failedDist = 0;
minDistMatrice = ORG_STRUC.minDistMatrice;
tic

% The variable to store number of structures generated with wrong space group:
ORG_STRUC.wrong_spacegroups = 0;

while goodPop < ORG_STRUC.populationSize + 0.5
    failedTime = toc;
    if (failedDist > 10000) || (failedTime > 300)
        if minDistMatrice(1,1) > 0.8*ORG_STRUC.minDistMatrice(1,1)
            if failedTime > 300
                %disp('WARNING! Can not generate a structure after 5 minutes. The minimum distance threshold will be lowered by 0.1.');
                USPEXmessage(504,'',0);
            else
                %disp('WARNING! Can not generate a structure after 10000 tries. The minimum distance threshold will be lowered by 0.1.');
                USPEXmessage(505,'',0);
            end
            
            minDistMatrice = 0.9*minDistMatrice;
            failedDist = 0;
            tic
        else
            disp('Could not generate a structure after 30000 tries or 15 minutes.');
            disp('Please check the input files. The calculation has to stop.');
            disp('Possible reasons:  unreasonably big contraints (IonDistances).');
            disp('Remember they should be much smaller than the real interatomic distances');
            disp('but not too small for pseudopotential overlap errors to kill interatomic repulsion.');
            quit;
        end
    end
    
    errorS = 0;
    Stokes = 0;
    candidate = rand(sum(ORG_STRUC.numIons),3);
    numIons = ORG_STRUC.numIons;
    POP_STRUC.POPULATION(goodPop).numIons = numIons;
    
    if ORG_STRUC.nsymN(1,1) == 0  % H. Stokes code to create a crystal with given symmetry
        sym_coef = ORG_STRUC.sym_coef;
        Stokes = 1;
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
        
        if ORG_STRUC.constLattice    % lat1 = lattice
            lat1 = ORG_STRUC.lattice;
        else                         % lat1 = lattice volume
            lat1 = ORG_STRUC.latVolume;
        end
        
        if sum(ORG_STRUC.splitInto) < 4
            cd(['CalcFoldTemp']);
            [candidate, lat, errorS] = symope_crystal(nsym, numIons, lat1, minDistMatrice, ORG_STRUC.sym_coef);
            
            Ind_No = goodPop;
            if errorS == 0
                [resulted_sg_no, nothing] = determine_spacegroup(lat, candidate, ORG_STRUC.atomType, ...
                    numIons, nsym, POP_STRUC.bodyCount + Ind_No, ...
                    ORG_STRUC.SGtolerance, ORG_STRUC.homePath);
            end
            
            cd(ORG_STRUC.homePath);
        else
            startLat = rand(6,1);
            startLat(4:6) = (pi/2);
            check_startLat = latConverter(startLat);
            volLat = det(check_startLat);
            ratio = ORG_STRUC.latVolume/volLat;
            startLat(1:3) = startLat(1:3)*(ratio)^(1/3);
            [lat, errorS, candidate] = splitBigCell(startLat, ORG_STRUC.splitInto(ceil(length(ORG_STRUC.splitInto)*rand)),numIons, nsym);
        end
    end
    
    if errorS == 1
        goodBad = 0;
    else
        goodBad = distanceCheck(candidate, lat, numIons, minDistMatrice*sym_coef);
    end
    
    if goodBad
        if nsym >= 1
            if nsym ~= resulted_sg_no
                ORG_STRUC.wrong_spacegroups = ORG_STRUC.wrong_spacegroups + 1;
            end
        end
        
        POP_STRUC.POPULATION(goodPop).LATTICE = lat;
        POP_STRUC.POPULATION(goodPop).COORDINATES = candidate;
        POP_STRUC.POPULATION(goodPop).howCome = '  Random  ';
        
        tic
        failedDist = 0;
        newSym = 1;
        
        if Stokes == 1
            disp(['Crystal ' num2str(goodPop) ' built with the symmetry group ' num2str(nsym) ' (' spaceGroups(nsym) ')'])
        else
            disp(['Structure ' num2str(goodPop) ' generated successfully'])
        end
        
        goodPop = goodPop + 1;
    else
        failedDist = failedDist + 1;
    end
end

if ORG_STRUC.wrong_spacegroups > 0
    disp([' ']);
    disp(['ATTENTION! In ' num2str(ORG_STRUC.wrong_spacegroups) ' / ' ...
        num2str(ORG_STRUC.initialPopSize) ' cases actually generated symmetry was different.']);
    disp([' ']);
end

%%%%%%%%%%%% END defining the first generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% SEEDING %%%%%%%%%%%%%%%%%
pick_Seeds();
if ORG_STRUC.doFing
    pickAntiSeeds();
end
Start_POP_PSO();
