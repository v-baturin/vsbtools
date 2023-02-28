function [numIons, offspring,potentialLattice,fracFrac,dimension,offset,fracLattice] = heredity_GB(par_one,par_two)

% The heredity was designed to allow structural heritage. 
% This function looks terribly complicated - but don't worry, it's all
% quite straight forward.
% Keep in mind that we are dealing with a periodical problem (periodicity : 1)
% => example (coordinate) : 2.41 = 1.41 = 0.41 = -0.59 = - 1.59...

global POP_STRUC
global ORG_STRUC
ordering_on = ORG_STRUC.ordering;

 lattice = POP_STRUC.POPULATION(par_one).GB_LATTICE;
 parent1 = POP_STRUC.POPULATION(par_one).GB_COORDINATES;
 par_order1 = POP_STRUC.POPULATION(par_one).GB_order;
 par_numIons1 = POP_STRUC.POPULATION(par_one).GB_numIons;
 atype1 = POP_STRUC.POPULATION(par_one).GB_typesAList;

 parent2 = POP_STRUC.POPULATION(par_two).GB_COORDINATES;
 par_order2 = POP_STRUC.POPULATION(par_two).GB_order;
 par_numIons2 = POP_STRUC.POPULATION(par_two).GB_numIons;
 atype2 = POP_STRUC.POPULATION(par_two).GB_typesAList;

 potentialLattice=lattice;
 numIons = ORG_STRUC.numIons; 
% fracFrac determines what spatial fraction of the one parent will be taken.
% The rest is taken from the other parent. Thus: 0.25 means one spatial
% fourth versus three fourths.

% fracFrac = 0.5;
fracFrac = 0.25 +rand(1)/2;
% in this way the most extreme difference is 1:3
% In the following lines we simply change from colomn vectors to matrices,
% where each colomn represents the coordinates on one lattice vector
% (it makes life easier)
%parent1 = POP_STRUC.POPULATION(par_one).COORDINATES;
%parent2 = POP_STRUC.POPULATION(par_two).COORDINATES;
% We select randomly a dimension used for the spacial criteria of the
% heredity
      dimension = RandInt(1,1,[1,2]);

 if rand(1)>ORG_STRUC.percSliceShift

  offset = rand(2,1);
  parent1(:,dimension)= parent1(:,dimension)+offset(1,1);
  parent2(:,dimension)= parent2(:,dimension)+offset(2,1);
  parent1(:,dimension) = parent1(:,dimension) - floor(parent1(:,dimension));
  parent2(:,dimension) = parent2(:,dimension) - floor(parent2(:,dimension));

 else
% The offset allows
  offset = rand(6,1);
  parent1(:,1)= parent1(:,1)+offset(1,1);
  parent1(:,2)= parent1(:,2)+offset(2,1);
  parent1(:,3)= parent1(:,3)+offset(3,1);
  parent2(:,1)= parent2(:,1)+offset(4,1);
  parent2(:,2)= parent2(:,2)+offset(5,1);
  parent2(:,3)= parent2(:,3)+offset(6,1);
% the offset will chuck some ions out of the [0,1] boundaries. We need to
% adjust this (otherwise the offset makes no sense at all)
% IMPORTANT : translation of all ions of one structure (parent) results in
% an equivalent structure. The reason for this translation is,
% that the localization of substructures change by inheritance (inhibiting convergence)
  parent1 = parent1 - floor(parent1);
  parent2 = parent2 - floor(parent2);
 end

% find the indices of all ions located in the right respective fraction of
% each structure
 whichIons1 = find(parent1(:,dimension)<fracFrac);
 whichIons2 = find(parent2(:,dimension)>fracFrac);

% chooses the more ordered slab out of 2 for both parents
 if ordering_on
  if length(whichIons1) > 0
   ord1 = sum(par_order1(whichIons1))/length(whichIons1);
  else
   ord1 = 0;
  end
  if length(whichIons2) > 0
   ord2 = sum(par_order2(whichIons2))/length(whichIons2);
  else
   ord2 = 0;
  end
  for i = 1 : 3  %Nslabs=3 for some reason
   offset_extra = rand(2,1);
   parent1(:,dimension) = parent1(:,dimension) + offset_extra(1,1);
   parent2(:,dimension) = parent2(:,dimension) + offset_extra(2,1);
   parent1(:,dimension) = parent1(:,dimension) - floor(parent1(:,dimension));
   parent2(:,dimension) = parent2(:,dimension) - floor(parent2(:,dimension));
   whichIons1_extra = find(parent1(:,dimension)<fracFrac);
   whichIons2_extra = find(parent2(:,dimension)>fracFrac);
   if length(whichIons1_extra) > 0
    ord1_extra = sum(par_order1(whichIons1_extra))/length(whichIons1_extra);
   else
    ord1_extra = 0;
   end
   if length(whichIons2_extra) > 0
    ord2_extra = sum(par_order2(whichIons2_extra))/length(whichIons2_extra);
   else
    ord2_extra = 0;
   end
   if ord1 > ord1_extra
    parent1(:,dimension) = parent1(:,dimension) - offset_extra(1,1);
    parent1(:,dimension) = parent1(:,dimension) - floor(parent1(:,dimension));
   else
    whichIons1 = whichIons1_extra;
    ord1 = ord1_extra;
   end
   if ord2 > ord2_extra
    parent2(:,dimension) = parent2(:,dimension) - offset_extra(2,1);
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

 ionCount = zeros(length(ORG_STRUC.atomType),1);

for ind=1:length(par_numIons1)
     ionChange1(ind)=sum(par_numIons1(1:ind));
end
ionCh1= zeros (1,length(ionChange1)+1);
ionCh1(2:end) = ionChange1;

for ind = 1:length(par_numIons2)
    ionChange2(ind) = sum(par_numIons2(1:ind));
end
ionCh2 = zeros (1,length(ionChange2)+1);
ionCh2(2:end) = ionChange2;

 for ind = 1:length(ORG_STRUC.atomType)
     if ORG_STRUC.numIons(ind)~=0
     % this loop is for all the different ion types.
     candidates = [];
     % find the indices of the ions for the type of ion done in this turn of
     % the loop
     smaller1 = find(whichIons1<=ionCh1(ind+1));
     smaller2 = find(whichIons2<=ionCh2(ind+1));
     ions1  = find(whichIons1(smaller1)>ionCh1(ind));
     ions2  = find(whichIons2(smaller2)>ionCh2(ind));
    % total number of this type of ion
     ionCount(ind,1) = length(ions1) + length(ions2);
     numatom = ORG_STRUC.numIons(ind);
	if ionCount(ind,1) > numatom
          howmany = ionCount(ind,1) - numatom;
        % we determine the amount to be eliminated from each fraction
        % probabilisticly (probability of each ion to be of a certain 
        %fraction linear with the size of this fraction)
        deleteIons_1 = length(find(rand(howmany,1)<fracFrac));
        deleteIons_2 = howmany - deleteIons_1;
        if length(ions1) < deleteIons_1
            oops = length(ions1) - deleteIons_1;
            deleteIons_1 = length(ions1);
            deleteIons_2 = deleteIons_2-oops;
        elseif length(ions2) < deleteIons_2
            oops = length(ions2) - deleteIons_2;
            deleteIons_2 = length(ions2);
            deleteIons_1 = deleteIons_1-oops;
        end
        if ordering_on
	    if isempty(ions1)
               ions1 = [];
            else
               atoms = sort_order_surface(par_order1, ions1, 'ascend');
               ions1(atoms(1:deleteIons_1)) = [];
            end
      	    if isempty(ions2)
               ions2 = [];
            else
	       atoms = sort_order_surface(par_order2, ions2, 'ascend');
               ions2(atoms(1:deleteIons_2)) = [];
            end
        else 
        % now we delete the surplus ions (randomly within each fraction) 
          for xy = 1:deleteIons_1
            ions1(RandInt(1,1,[1,length(ions1)])) = [];
          end
          for xy = 1:deleteIons_2
            ions2(RandInt(1,1,[1,length(ions2)])) = [];
          end
        end
            ionCount(ind,1)=numatom;
%            disp(['too many atoms: delete']);
    elseif ionCount(ind,1) < numatom
%%%%%%%%%%%%%%Add atoms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       disp(['not enough atoms: add']);
        howmany = numatom - ionCount(ind,1);
        % Linear probability of each ion to be of a fraction (linear with the size of the fraction)
        suppleIons_1 = length(find(rand(howmany,1)<fracFrac));
        suppleIons_2 = howmany-suppleIons_1;
        tempInd = find(parent2(ionCh2(ind)+1:ionCh2(ind+1),dimension)<fracFrac);
           if length(tempInd)==0
               suppleIons_1 = 0;
               suppleIons_2 = howmany;
           end
        % for each fraction we choose the determined amount of ions. The
        % new ions are chosen (randomly) from the structure the ions in this fraction
        % weren't originally chosen from
        if suppleIons_1 > 0
           tempInd = find(parent2(ionCh2(ind)+1:ionCh2(ind+1),dimension)<fracFrac);
           if isempty(tempInd)
               disp('');
               parent2
           end
           if ordering_on
            atoms = sort_order_surface(par_order2, tempInd, 'descend');
            candidates = parent2(ionCh2(ind) + tempInd(atoms(1:suppleIons_1)) , :);
           else
            crutch = rand(length(tempInd),1);
            [junk, newOrder] = sort(crutch);
            candidates = parent2(ionCh2(ind) + tempInd(newOrder(1:suppleIons_1)) , :);
           end
        end
        if suppleIons_2 > 0
            tempInd = find(parent1(ionCh1(ind)+1:ionCh1(ind+1),dimension)>fracFrac);
            if isempty(tempInd)
               disp('');
               parent1
            end

           if ordering_on
            atoms = sort_order_surface(par_order1, tempInd, 'descend');
            candidates = cat(1,candidates,parent1(ionCh1(ind) + tempInd(atoms(1:suppleIons_2)) , :));
           else
            crutch = rand(length(tempInd),1);
            [junk, newOrder] = sort(crutch);
            candidates = cat(1,candidates,parent1(ionCh1(ind)+ tempInd(newOrder(1:suppleIons_2)) , :));
           end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fill the selected ions coordinates into the offspring
    addOn = cat(1,parent1(whichIons1(smaller1(ions1(:))) , :), parent2(whichIons2(smaller2(ions2(:))) , :));
    addOn = cat(1,addOn,candidates);
   
    offspring = cat(1,offspring,addOn);
     end
 end
