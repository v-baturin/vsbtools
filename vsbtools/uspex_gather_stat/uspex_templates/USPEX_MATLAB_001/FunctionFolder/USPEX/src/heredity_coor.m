function offspring = heredity_coor(parent1, whichIons1, order1, numIons1,...
                                   parent2, whichIons2, order2, numIons2,...
                                 numIons, frac, dimension, atomType, cor_dir)

% find the indices of all ions located in the right respective
ionCh1 = heredity_ionCh(numIons1);
ionCh2 = heredity_ionCh(numIons2);

ionCount = zeros(length(atomType),1);
offspring = [];
for ind = 1 : length(atomType)  % loop through all different ion types

    candidates = [];
    % find the indices of the ions of the type ind
    smaller1 = find(whichIons1<=ionCh1(ind+1));
    smaller2 = find(whichIons2<=ionCh2(ind+1));
    ions1 = find(whichIons1(smaller1)>ionCh1(ind));
    ions2 = find(whichIons2(smaller2)>ionCh2(ind));

    % total number of this type of ion
    ionCount(ind,1) = length(ions1) + length(ions2);

    % Add if the number of ions is too small
    if ionCount(ind,1) < numIons(ind)

        Addhowmany   = numIons(ind) - ionCount(ind,1);
        suppleIons_1 = length(find(rand(Addhowmany,1)<frac));
        suppleIons_2 = Addhowmany - suppleIons_1;
        if suppleIons_1 + length(ions2) > numIons2(ind)
           suppleIons_1 = numIons2(ind) - length(ions2);
           suppleIons_2 = Addhowmany - suppleIons_1;
        elseif suppleIons_2 + length(ions1) > numIons1(ind)
           suppleIons_2 = numIons1(ind) - length(ions1);
           suppleIons_1 = Addhowmany - suppleIons_2;
        end
        tmp =  heredity_Add(suppleIons_1, parent2, order2, ionCh2, ind,...
                                              cor_dir, frac, dimension, 1);
        candidates = [candidates; tmp];
        tmp =  heredity_Add(suppleIons_2, parent1, order1, ionCh1, ind,...
                                              cor_dir, frac, dimension, 2);
        candidates = [candidates; tmp];

    % Delete if the number of ions is too big
    elseif ionCount(ind,1) > numIons(ind)

        Delhowmany   = ionCount(ind,1) - numIons(ind);
        deleteIons_1 = length(find(rand(Delhowmany,1)<frac));
        deleteIons_2 = Delhowmany - deleteIons_1;
        % quick safeguarding against trouble
        if length(ions1) < deleteIons_1
            oops = deleteIons_1 - length(ions1);
            deleteIons_1 = length(ions1);
            deleteIons_2 = deleteIons_2 + oops;
        elseif length(ions2) < deleteIons_2
            oops = deleteIons_2 - length(ions2);
            deleteIons_2 = length(ions2);
            deleteIons_1 = deleteIons_1 + oops;
        end

        ions1 = heredity_Delete(ions1, order1, deleteIons_1, cor_dir);
        ions2 = heredity_Delete(ions2, order2, deleteIons_2, cor_dir);
    end

    % fill the selected ions coordinates into the offspring
    offspring = [offspring; parent1(whichIons1(smaller1(ions1(:))) , :)];
    offspring = [offspring; parent2(whichIons2(smaller2(ions2(:))) , :)];
    offspring = [offspring; candidates];
end

%-------------------------------------------------------------------------
function ionCh = heredity_ionCh(numIons)
%-------------------------------------------------------------------------
for ind = 1 : length(numIons)
    ionChange(ind) = sum(numIons(1:ind));
end
ionCh = zeros(1,length(ionChange)+1);
ionCh(2:end) = ionChange;
                                                                                
%-------------------------------------------------------------------------
function ions = heredity_Delete(ions, order, delete_Ions, cor_dir)
%-------------------------------------------------------------------------
% Delete the surplus ions (either randomly or the least ordered ones) 
if length(ions) > 0
   if cor_dir <= 0  % low fitness (=good) <=> high order
      atoms = sort_order(order, ions, 'ascend');
   else
      atoms = sort_order(order, ions, 'descend');
   end
   ions(atoms(1:delete_Ions)) = [];
end

%-------------------------------------------------------------------------
function ToAdd = heredity_Add(num, parent, order, ionCh, ...
                           ind, cor_dir, frac, dimension, ID)
%-------------------------------------------------------------------------
%disp(['Adding: ' num2str([num, ind, ionCh(ind)+1, ionCh(ind+1)]) ]);
if num > 0
   if ID == 1
       tempInd = find(parent(ionCh(ind)+1:ionCh(ind+1),dimension)<frac);
   else
       tempInd = find(parent(ionCh(ind)+1:ionCh(ind+1),dimension)>frac);
   end

   if cor_dir <= 0  % low fitness (=good) <=> high order
      atoms = sort_order(order, tempInd, 'descend');
   else % low fitness (=good) <=> low order
      atoms = sort_order(order, tempInd, 'ascend');
   end
   tmp   = tempInd(atoms(1:num));
   ToAdd = parent(ionCh(ind) + tmp , :);
else
   ToAdd = [];
end

