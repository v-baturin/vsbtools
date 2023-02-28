function [numMols, numBlocks, offspring, potentialLattice, fracFrac,dimension,offset,MolTable] = heredity_molecule_var(par_one,par_two)

% The heredity was designed to allow structural heritage. This is achieved
% by

% This function looks terribly complicated - but don't worry, it's all
% quite straight forward.
% Keep in mind that we are dealing with a periodical problem (periodicity : 1)
% => example (coordinate) : 2.41 = 1.41 = 0.41 = -0.59 = - 1.59...

global POP_STRUC
global ORG_STRUC

ordering_on = ORG_STRUC.ordering;

% fracFrac determines what spatial fraction of the one parent will be taken.
% The rest is taken from the other parent. Thus: 0.25 means one spatial
% fourth versus three fourths.
MolTable = [];
fracFrac = 0.25 + rand(1)/2;
% in this way the most extreme difference is 1:3
numIons1 = POP_STRUC.POPULATION(par_one).numMols;
numIons2 = POP_STRUC.POPULATION(par_two).numMols;

parent1 = zeros(sum(numIons1),3);
parent2 = zeros(sum(numIons2),3);

for ind = 1 : sum(numIons1)
    parent1(ind,:) = POP_STRUC.POPULATION(par_one).MOLECULES(ind).MOLCENTER;
    order1(ind)    = POP_STRUC.POPULATION(par_one).MOLECULES(ind).order;
end

for ind = 1 : sum(numIons2)
    parent2(ind,:) = POP_STRUC.POPULATION(par_two).MOLECULES(ind).MOLCENTER;
    order2(ind)    = POP_STRUC.POPULATION(par_two).MOLECULES(ind).order;
end

lat1 = POP_STRUC.POPULATION(par_one).LATTICE;
lat2 = POP_STRUC.POPULATION(par_two).LATTICE;
parent1 = parent1/lat1;
parent2 = parent2/lat2;

dimension = RandInt(1,1,[1,3]);

correlation_coefficient = ORG_STRUC.correlation_coefficient;
cor_dir = ORG_STRUC.cor_dir;

% determine the thickness of the unit cell while cutting in a given direction (not equal to the lattice vector for non-cubic cells!)
if dimension == 1 
 L1 = abs(det(lat1)/norm(cross(lat1(2,:),lat1(3,:))));
 L2 = abs(det(lat2)/norm(cross(lat2(2,:),lat2(3,:))));
elseif dimension == 2
 L1 = abs(det(lat1)/norm(cross(lat1(1,:),lat1(3,:))));
 L2 = abs(det(lat2)/norm(cross(lat2(1,:),lat2(3,:))));
else
 L1 = abs(det(lat1)/norm(cross(lat1(1,:),lat1(2,:))));
 L2 = abs(det(lat2)/norm(cross(lat2(1,:),lat2(2,:))));
end
Lchar1 = 0.5*power(abs(det(lat1))/sum(numIons1), 1/3); % characteristic length - approximate 'radius' of the molecule in the cell = 0.5*(V/N)^1/3
Lchar2 = 0.5*power(abs(det(lat2))/sum(numIons2), 1/3);

Nslabs1 = round(L1/(Lchar1+(L1-Lchar1)*(cos(correlation_coefficient*pi/2))^2));
Nslabs2 = round(L2/(Lchar2+(L2-Lchar2)*(cos(correlation_coefficient*pi/2))^2));


if rand(1) > ORG_STRUC.percSliceShift    % shift in one dimension
 offset = rand(2,1);
 parent1(:,dimension) = parent1(:,dimension) + offset(1,1);
 parent2(:,dimension) = parent2(:,dimension) + offset(2,1);
 parent1(:,dimension) = parent1(:,dimension) - floor(parent1(:,dimension));
 parent2(:,dimension) = parent2(:,dimension) - floor(parent2(:,dimension));

else     % shift in all dimensions

offset = rand(6,1);
parent1(:,1)= parent1(:,1)+offset(1,1);
parent1(:,2)= parent1(:,2)+offset(2,1);
parent1(:,3)= parent1(:,3)+offset(3,1);
parent2(:,1)= parent2(:,1)+offset(4,1);
parent2(:,2)= parent2(:,2)+offset(5,1);
parent2(:,3)= parent2(:,3)+offset(6,1);
parent1 = parent1 - floor(parent1);
parent2 = parent2 - floor(parent2);

end

% find the indices of all ions located in the right respective fraction of each structure
whichIons1 = find(parent1(:,dimension)<fracFrac);
whichIons2 = find(parent2(:,dimension)>fracFrac);

% chooses the more ordered slab out of 2 for both parents%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ordering_on
 if length(whichIons1) > 0
  ord1 = sum(order1(whichIons1))/length(whichIons1);
 elseif cor_dir <= 0
  ord1 = 0;
 else
  ord1 = 1;
 end
 if length(whichIons2) > 0
  ord2 = sum(order2(whichIons2))/length(whichIons2);
 elseif cor_dir <= 0
  ord2 = 0;
 else
  ord2 = 1;
 end
 for i = 1 : Nslabs1
  offset_extra = rand;
  parent1(:,dimension) = parent1(:,dimension) + offset_extra;
  parent1(:,dimension) = parent1(:,dimension) - floor(parent1(:,dimension));
  whichIons1_extra = find(parent1(:,dimension)<fracFrac);
  if length(whichIons1_extra) > 0
   ord1_extra = sum(order1(whichIons1_extra))/length(whichIons1_extra);
  elseif cor_dir <= 0
   ord1_extra = 0;
  else
   ord1_extra = 1;
  end
  if ((ord1 > ord1_extra) & (cor_dir <= 0)) | ((ord1 < ord1_extra) & (cor_dir == 1))
   parent1(:,dimension) = parent1(:,dimension) - offset_extra;
   parent1(:,dimension) = parent1(:,dimension) - floor(parent1(:,dimension));
  else
   whichIons1 = whichIons1_extra;
   ord1 = ord1_extra;
  end
 end
 for i = 1 : Nslabs2
  offset_extra = rand;
  parent2(:,dimension) = parent2(:,dimension) + offset_extra;
  parent2(:,dimension) = parent2(:,dimension) - floor(parent2(:,dimension));
  whichIons2_extra = find(parent2(:,dimension)>fracFrac);
  if length(whichIons2_extra) > 0
   ord2_extra = sum(order2(whichIons2_extra))/length(whichIons2_extra);
  elseif cor_dir <= 0 
   ord2_extra = 0;
  else
   ord2_extra = 1;
  end
  if ((ord2 > ord2_extra) & (cor_dir <= 0)) | ((ord2 < ord2_extra) & (cor_dir == 1))
   parent2(:,dimension) = parent2(:,dimension) - offset_extra;
   parent2(:,dimension) = parent2(:,dimension) - floor(parent2(:,dimension));
  else
   whichIons2 = whichIons2_extra;
   ord2 = ord2_extra;
  end
 end
end

% prepare the variable offspring. This needs to be done due to the command
% 'cat' used later on
offspring = zeros(0,3);
ionCount = zeros(size(ORG_STRUC.numMols,2),1);

  for ind = 1 : length(numIons1)
    ionChange1(ind) = sum(numIons1(1:ind));
  end
  ionCh1 = zeros(1,length(ionChange1)+1);
  ionCh1(2:end) = ionChange1;
  for ind = 1 : length(numIons2)
    ionChange2(ind) = sum(numIons2(1:ind));
  end
  ionCh2 = zeros(1,length(ionChange2)+1);
  ionCh2(2:end) = ionChange2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we have to fix the correct number of atoms/blocks in the offspring %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 maxBlocks = POP_STRUC.POPULATION(par_one).numMols/ORG_STRUC.numMols + POP_STRUC.POPULATION(par_two).numMols/ORG_STRUC.numMols;
 child = zeros(1, size(ORG_STRUC.numIons,2));
 for i = 1 : size(ORG_STRUC.numIons,2)
    smaller1 = find(whichIons1<=ionCh1(i+1));
    smaller2 = find(whichIons2<=ionCh2(i+1));
    ions1 = find(whichIons1(smaller1)>ionCh1(i));
    ions2 = find(whichIons2(smaller2)>ionCh2(i));
    child(i) = length(ions1) + length(ions2);
 end
 [desired_comp, desired_blocks] = findDesiredComposition(maxBlocks, ORG_STRUC.numIons, child);

for ind = 1 : size(ORG_STRUC.numMols,2)  % loop through all different ion types

    candidates = []; %to store the supplementary coordinates
    supp = [];       %to store the supplementary tables
    % find the indices of the ions of the type ind
    smaller1 = find(whichIons1<=ionCh1(ind+1));
    smaller2 = find(whichIons2<=ionCh2(ind+1));
    ions1 = find(whichIons1(smaller1)>ionCh1(ind));
    ions2 = find(whichIons2(smaller2)>ionCh2(ind));

    % total number of this type of ion
    ionCount(ind,1) = length(ions1) + length(ions2);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if the number of ions we have is smaller than the number we need, then we have to fill up with further ions

    if ionCount(ind,1) < desired_comp(ind)

        % calculate how many ions (of this type) are missing
        howmany = desired_comp(ind) - ionCount(ind,1);

        % determine how many ions are to be filled into each fraction probabilistically.
        % Linear probability of each ion to be of a fraction (linear with the size of the fraction)
        % varcomp: make sure we have enough atoms in each parent stucture, if not - adjust suppleIon


        suppleIons_1 = length(find(rand(howmany,1)<fracFrac));
        suppleIons_2 = howmany - suppleIons_1;
          temp2 = length(find(parent2(ionCh2(ind)+1:ionCh2(ind+1),dimension)<fracFrac));
          temp1 = length(find(parent1(ionCh1(ind)+1:ionCh1(ind+1),dimension)>fracFrac));
          if suppleIons_1 > temp2
            suppleIons_2 = suppleIons_2 + (suppleIons_1 - temp2);
            suppleIons_1 = temp2;
          end
          if suppleIons_2 > temp1
            suppleIons_1 = suppleIons_1 + (suppleIons_2 - temp1);
            suppleIons_2 = temp1;
          end

        % for each fraction we choose the determined amount of ions. The new ions are chosen (randomly) 
        % from the structure the ions in this fraction weren't originally chosen from
           supp=zeros(suppleIons_1+suppleIons_2,3);
        if suppleIons_1 > 0
           tempInd = find(parent2(ionCh2(ind)+1:ionCh2(ind+1),dimension)<fracFrac);
           if ordering_on
            if cor_dir <= 0  % low fitness (=good) <=> high order
              atoms = sort_order(order2, tempInd, 'descend');
            else % low fitness (=good) <=> low order
              atoms = sort_order(order2, tempInd, 'ascend');
            end
            candidates = parent2(ionCh2(ind) + tempInd(atoms(1:suppleIons_1)) , :);
            supp(1:suppleIons_1,2)=ionCh2(ind) + tempInd(atoms(1:suppleIons_1));

           else
            crutch = rand(length(tempInd),1);
            [junk, newOrder] = sort(crutch);
            candidates = parent2(ionCh2(ind) + tempInd(newOrder(1:suppleIons_1)) , :);
            supp(1:suppleIons_1,2)=ionCh2(ind) + tempInd(newOrder(1:suppleIons_1));
           end
            supp(1:suppleIons_1,1)=2;
            supp(1:suppleIons_1,3)=ind;
        end

        if suppleIons_2 > 0
           tempInd = find(parent1(ionCh1(ind)+1:ionCh1(ind+1),dimension)>fracFrac);
           if ordering_on
            if cor_dir <= 0  % low fitness (=good) <=> high order
              atoms = sort_order(order1, tempInd, 'descend');
            else % low fitness (=good) <=> low order
              atoms = sort_order(order1, tempInd, 'ascend');
            end
            candidates = cat(1,candidates,parent1(ionCh1(ind) + tempInd(atoms(1:suppleIons_2)) , :));
            supp(suppleIons_1+1:suppleIons_1+suppleIons_2,2)=ionCh1(ind) + tempInd(atoms(1:suppleIons_2));
           else
            crutch = rand(length(tempInd),1);
            [junk, newOrder] = sort(crutch);
            candidates = cat(1,candidates,parent1(ionCh1(ind)+ tempInd(newOrder(1:suppleIons_2)) , :));
            supp(suppleIons_1+1:suppleIons_1+suppleIons_2,2)=ionCh1(ind) + tempInd(newOrder(1:suppleIons_2));
           end
           supp(suppleIons_1+1:suppleIons_1+suppleIons_2,1)=1;
           supp(suppleIons_1+1:suppleIons_1+suppleIons_2,3)=ind;
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if the number of ions we have is too big, we need to eliminate a few %

    elseif ionCount(ind,1) > desired_comp(ind)

        howmany = ionCount(ind,1) - desired_comp(ind);

        % we determine the amount to be eliminated from each fraction probabilisticly 
        % probability of each ion to be of a certain is fraction linear with the size of this fraction

        deleteIons_1 = length(find(rand(howmany,1)<fracFrac));
        deleteIons_2 = howmany - deleteIons_1;

        % quick safeguarding against trouble
        if length(ions1) < deleteIons_1
            oops = length(ions1) - deleteIons_1;
            deleteIons_1 = length(ions1);
            deleteIons_2 = deleteIons_2 - oops;
        elseif length(ions2) < deleteIons_2
            oops = length(ions2) - deleteIons_2;
            deleteIons_2 = length(ions2);
            deleteIons_1 = deleteIons_1 - oops;
        end

        % now we delete the surplus ions (either randomly or the least ordered ones) 
        if ordering_on
         if length(ions1) > 0
          if cor_dir <= 0  % low fitness (=good) <=> high order
            atoms = sort_order(order1, ions1, 'ascend');
          else
            atoms = sort_order(order1, ions1, 'descend');
          end
          ions1(atoms(1:deleteIons_1)) = [];
         end
         if length(ions2) > 0
          if cor_dir <= 0  % low fitness (=good) <=> high order
            atoms = sort_order(order2, ions2, 'ascend');
          else
            atoms = sort_order(order2, ions2, 'descend');
          end
          ions2(atoms(1:deleteIons_2)) = [];
         end
        else
          for xy = 1 : deleteIons_1
            ions1(RandInt(1,1,[1,length(ions1)])) = [];
          end
          for xy = 1 : deleteIons_2
            ions2(RandInt(1,1,[1,length(ions2)])) = [];
          end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fill the selected ions coordinates into the offspring
    addOn = cat(1,parent1(whichIons1(smaller1(ions1(:))) , :), parent2(whichIons2(smaller2(ions2(:))) , :));
    addOn = cat(1,addOn,candidates);

    offspring = cat(1,offspring,addOn);
    table1=[];
       if ~isempty(smaller1(ions1(:)))
          for i=1:length(smaller1(ions1(:)))
              table1(i,1) = 1;  %Parent Identity
              table1(i,2) = whichIons1(smaller1(ions1(i))); %Molecule index
              table1(i,3) = ind;  %Molecule type
          end
       end
       if ~isempty(smaller2(ions2(:)))
          for i=1:length(smaller2(ions2(:)))
              table1(length(smaller1(ions1(:)))+i,1)=2;
              table1(length(smaller1(ions1(:)))+i,2)=whichIons2(smaller2(ions2(i)));
              table1(length(smaller1(ions1(:)))+i,3)=ind;
          end
       end
%       if ~isempty(supp)
          table1=cat(1,table1,supp);
%       end
       MolTable=cat(1,MolTable,table1); 
   end

  numMols = desired_comp;
  numBlocks = desired_blocks;


%%%%%%%%%%%%%%%%%%%%%%%%%%
% % substract the offset (thus returning to the true values)
% offspring(:,dimension) = offspring(:,dimension) - offset;
% % adjust to the [0,1] boundaries
% offspring(:,dimension) = offspring(:,dimension)  - floor(offspring(:,dimension) );
%
% HOWEVER: substracting the offset is not necessary, since it changes
% nothing regarding the structure. If we leave the offset as it is, then we
% leave the substructures at their changed position. I asume this to be
% good against convergence.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~ORG_STRUC.constLattice

    fracLattice = rand(1);
    temp_potLat = fracLattice*lat1 + (1-fracLattice)*lat2;
    volLat = det(temp_potLat);
    if sign(volLat)==-1
        temp_potLat = -1*temp_potLat;
    end

    % scale the lattice to the volume we asume it approximately to be
    if (size(ORG_STRUC.numIons,1) == 2)
       N_T = size(ORG_STRUC.numMols,2);
       ch = POP_STRUC.convex_hull;
       try
     % works only for 2 types of atoms!
         x = numMols(N_T)/sum(numMols(1:N_T));
         i = 1;
         while x > ch([i+1],N_T)/sum(ch([i+1],1:N_T))
          i = i + 1;
         end
         v1 = det(POP_STRUC.POPULATION(ch([i+1],N_T+2)).LATTICE);
         v2 = det(POP_STRUC.POPULATION(ch([i],N_T+2)).LATTICE);
         a1 = ch([i+1],N_T)/sum(ch([i+1],1:N_T));
         a2 = ch([i],N_T)/sum(ch([i],1:N_T));
         latVol = abs((x-a2)/(a2-a1))*v1 + abs((x-a1)/(a2-a1))*v2;
       catch
        latVol = 0;
        for it = 1 : length(ORG_STRUC.latVolume)
           latVol = latVol + ionCount(it)*ORG_STRUC.latVolume(it);
        end
       end
    else
        latVol = det(lat1)*fracFrac + det(lat2)*(1-fracFrac);
    end
    ratio = latVol/volLat;
    temp_potLat = latConverter(temp_potLat);
    temp_potLat(1:3)= temp_potLat(1:3)*(ratio)^(1/3);
    potentialLattice = latConverter(temp_potLat);
else
    potentialLattice=ORG_STRUC.lattice;
end
