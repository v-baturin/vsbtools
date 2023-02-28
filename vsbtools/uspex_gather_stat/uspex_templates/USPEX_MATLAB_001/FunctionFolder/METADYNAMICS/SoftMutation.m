function nextID = SoftMutation(OFF_ID, POP_ID, lat)
global POP_STRUC
global ORG_STRUC
global OFF_STRUC
%------------------------------------------------------
%-------Step1: Set up the inputs
%------------------------------------------------------
nextID  = OFF_ID;
coor    = POP_STRUC.POPULATION(POP_ID).COORDINATES;
numIons = POP_STRUC.POPULATION(POP_ID).numIons;
cell    = POP_STRUC.POPULATION(POP_ID).superCell;

lat = bsxfun(@times, lat, cell); %Matrix multiplies Vector
maxIncrease   = ORG_STRUC.maxIncrease; % maximum cell multiplication factor when supercell is built

%------------------------------------------------------
%------Step2: calculate which softmodes to start
%------------------------------------------------------
if POP_STRUC.POPULATION(POP_ID).Softmode_num == 0
   %New structure
   [freq, eigvector, supercells_sizes] = calcSoftModes_varcomp(...
                   lat, coor, numIons, ORG_STRUC.maxAt, maxIncrease);
   POP_STRUC.POPULATION(POP_ID).eignFre = freq;
   POP_STRUC.POPULATION(POP_ID).eignVec = eigvector;
   POP_STRUC.POPULATION(POP_ID).eignSupercells = supercells_sizes;
end
   freq             = POP_STRUC.POPULATION(POP_ID).eignFre;
   eigvector        = POP_STRUC.POPULATION(POP_ID).eignVec;
   supercells_sizes = POP_STRUC.POPULATION(POP_ID).eignSupercells;
   last_freq = POP_STRUC.POPULATION(POP_ID).Softmode_num;
   good_freq = ObtainFreq(lat, coor, numIons, freq, eigvector, supercells_sizes, last_freq);
%------------------------------------------------------
%------Step3: try soft mutation until it succeeds.
%------------------------------------------------------
for f = good_freq : length(freq)
    good1 = 0;
    good2 = 0;
    POP_STRUC.POPULATION(POP_ID).Softmode_num = f;
    [eigV, coord0, lat0, numIons0] = calcEigenvectorK(eigvector(:,f), supercells_sizes(f,:), lat, coor, numIons);
    superCellSize = supercells_sizes(f,:);
    % opposite direction is degenerate using fingerprints criteria
    opposite_degenerate = checkDegenerate(coord0, numIons0, lat0, eigV);
    %---SoftMutation: direction 1  
    [good1, MUT_LAT, MUT_COORD, deviation] = Do_soft_mutation(coord0, numIons0, lat0, eigV(:,1));
    if good1
        SaveOFF(POP_ID, OFF_ID, MUT_COORD, MUT_LAT, numIons0, freq, f, ...
                MUT_COORD-coord0, deviation, superCellSize);
        nextID = nextID + 1;
    else
       disp(['Direction 1 at mode ' num2str(f) ' is not good']);
       disp(['Mut_max = ' num2str(max(deviation))]);
    end
    %---SoftMutation: direction 2
    if opposite_degenerate == 0
       [good2, MUT_LAT, MUT_COORD, deviation] = Do_soft_mutation(coord0, numIons0, lat0, -1*eigV(:,1));
       if good2
          if good1 == 1
             ToWrite = length(OFF_STRUC.POPULATION)+1;
          else
             ToWrite = OFF_ID;
          end
          SaveOFF(POP_ID, ToWrite, MUT_COORD, MUT_LAT, numIons0, freq, f, ...
                  MUT_COORD-coord0, deviation, superCellSize);
          nextID = nextID + 1;
       else
           disp(['Direction 2 at mode ' num2str(f) ' is not good']);
           disp(['Mut_max = ' num2str(max(deviation))]);
       end
    end 

    %---------------------------------------------------------------------------------
    if max([good1, good2])>0
       break
    elseif f==length(freq)
       nextID = 0;  %We don't try from now, very IMPORTANT...
       disp('Out of soft modes.........OOPPPS..')
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   END creating mutants   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SaveOFF(POP_ID, OFF_ID, MUT_COORD, MUT_LAT, numIons0, freq, f, vec, deviation, superCellSize)
%Save the necessary items if the mutation is accepted
global POP_STRUC
global OFF_STRUC

OFF_STRUC.POPULATION(OFF_ID).COORDINATES = MUT_COORD;
OFF_STRUC.POPULATION(OFF_ID).LATTICE     = MUT_LAT;
OFF_STRUC.POPULATION(OFF_ID).numIons     = numIons0;
OFF_STRUC.POPULATION(OFF_ID).superCell   = superCellSize.*POP_STRUC.POPULATION(POP_ID).superCell;
info_parents = struct('parent', {},'mut_degree', {}, ...
     'mut_vec',{},'mut_mode',{},'mut_fre',{},'enthalpy',{});
info_parents(1).parent  = POP_STRUC.POPULATION(POP_ID).Number;
info_parents.enthalpy   = POP_STRUC.POPULATION(POP_ID).Enthalpies(end);       
info_parents.mut_degree = deviation;
info_parents.mut_mode   = f;
info_parents.mut_vec    = vec;
info_parents.mut_fre    = freq(f);
OFF_STRUC.POPULATION(OFF_ID).Parents = info_parents;
disp(['Structure ' num2str(OFF_ID) ' softmutated in freq = ' num2str(freq(f),'%6.3f') ...
      ' (mode # ' num2str(f) ')']);
disp(['Mut_max   = ' num2str(max(deviation),'%6.3f') ', Mut_average = ' num2str(mean(deviation),'%6.3f')]);
disp(['SuperCell = ' num2str(OFF_STRUC.POPULATION(OFF_ID).superCell)]); 
disp(['  ']);

function goodfreq = ObtainFreq(lat, coor, numIons, freq, eigvector, supercells_sizes, lastfreq);
%%%%%%%%%%% Check the degeneracy of frequencies using different criteria; 
%%%%%%%%%%% fingerprint criterion used only for non-stochastic softmutation
global POP_STRUC
global OFF_STRUC
global ORG_STRUC

weight        = ORG_STRUC.weight;
atomType      = ORG_STRUC.atomType;
toleranceFing = ORG_STRUC.toleranceFing;
if lastfreq > 0
   [eigV, coord0, lat0, numIons0]     = calcEigenvectorK(eigvector(:,lastfreq), ...
                                    supercells_sizes(lastfreq,:), lat, coor, numIons);
   [MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV, 1);
else %The structure itself --- to distinguish accoustic mode
    MUT_LAT0 = lat;
    MUT_COORD0 = coor;
    numIons0 = numIons;
end

[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
[order1, f1, atom_fing1]           = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);
goodfreq = lastfreq;
for i = lastfreq+1 : length(freq)
    goodfreq = goodfreq+1;
   if abs(freq(goodfreq)) > 1.0e-6   %we don't use zeros frequency
%%%%fing for the next one
      [eigV, coord0, lat0, numIons0]      = calcEigenvectorK(eigvector(:,i), supercells_sizes(i,:), lat, coor, numIons);
      [MUT_LAT0, MUT_COORD0, deviation0]  = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV, 1);
      [Ni, V, dist_matrix, typ_i, typ_j]  = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
      [order1, f2, atom_fing1]            = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);
%%%%fing for the last one
       dist = cosineDistance(f1, f2, weight);
    else
       dist = 0;
    end

    if (dist > toleranceFing)
       break
    else
       %disp('Find a degenerate mode, skip it.....')
    end
end

function degenerate = checkDegenerate(coord0, numIons0, lat0, eigV)
global ORG_STRUC
weight        = ORG_STRUC.weight;
atomType      = ORG_STRUC.atomType;
tolerance     = ORG_STRUC.toleranceFing;

[MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV, 1);
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
[order1, f1, atom_fing1]           = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);

[MUT_LAT0, MUT_COORD0, deviation0] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, -1*eigV, 1);
[Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(MUT_LAT0, MUT_COORD0, numIons0, atomType);
[order1, f2, atom_fing1]           = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons0);

dist = cosineDistance(f1, f2, weight);
if dist < tolerance
   degenerate = 1;
   disp('Opposite direction leads to same structure, skip it.....')
else
   degenerate = 0;
end

function  [goodAtomMutant, MUT_LAT, MUT_COORD, deviation] = Do_soft_mutation(coord0, numIons0, lat0, eigV)
global ORG_STRUC
for i = 0 : 10
   [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_Mutation(coord0, numIons0, lat0, eigV(:,1), 1-i/13);
   goodAtomMutant = distanceCheck(MUT_COORD, MUT_LAT, numIons0, ORG_STRUC.minDistMatrice);
   if goodAtomMutant 
      break
   end
end
