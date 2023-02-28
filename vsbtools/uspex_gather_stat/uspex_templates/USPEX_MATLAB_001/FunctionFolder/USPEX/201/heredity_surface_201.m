function [numIons, offspring,potentialLattice,fracFrac,dimension,offset] = heredity_surface_201(par_one,par_two,cell)

global POP_STRUC
global ORG_STRUC

ordering_on = ORG_STRUC.ordering;
atomType = ORG_STRUC.atomType;
lat_base = ORG_STRUC.bulk_lat;

lat1     = POP_STRUC.POPULATION(par_one).Surface_LATTICE;
coor1    = POP_STRUC.POPULATION(par_one).Surface_COORDINATES;
numIons1 = POP_STRUC.POPULATION(par_one).Surface_numIons;

lat2     = POP_STRUC.POPULATION(par_two).Surface_LATTICE;
coor2    = POP_STRUC.POPULATION(par_two).Surface_COORDINATES;
numIons2 = POP_STRUC.POPULATION(par_two).Surface_numIons;

[lat1, parent1, par_order1, par_numIons1]=cresupercell(lat1, coor1, numIons1, lat_base, cell, atomType);
[lat2, parent2, par_order2, par_numIons2]=cresupercell(lat2, coor2, numIons2, lat_base, cell, atomType);

% fracFrac determines what spatial fraction of the one parent will be taken.
% The rest is taken from the other parent. Thus: 0.25 means one spatial
% fourth versus three fourths.
numIons = [];
offspring = [];
potentialLattice = [];
fracFrac = [];
dimension = [];
offset = [];
if ~isempty(lat1) && ~isempty(lat2)
fracFrac = 0.25 +rand(1)/2;
if size(parent1,1)<2 && size(parent2,1)<2
   fracFrac=1;
   offset=0.0;
   dimension=3;
   offspring=zeros(0,3);
   offspring=cat(1,parent1,parent2); 
   numIons=par_numIons1+par_numIons2;
else
   dimension = RandInt(1,1,[1,2]);
end

correlation_coefficient = POP_STRUC.correlation_coefficient;
cor_dir = POP_STRUC.cor_dir;

% determine the thickness of the unit cell while cutting in a given direction (not equal to the lattice vector for non-cubic cells!)
if dimension == 1
    L1 = abs(det(lat1)/norm(cross(lat1(2,:),lat1(3,:))));
    L2 = abs(det(lat2)/norm(cross(lat2(2,:),lat2(3,:))));
else
    L1 = abs(det(lat1)/norm(cross(lat1(1,:),lat1(3,:))));
    L2 = abs(det(lat2)/norm(cross(lat2(1,:),lat2(3,:))));
end
Lchar1 = 0.5*power(abs(det(lat1))/sum(par_numIons1), 1/3); % characteristic length - approximate 'radius' of the atom in the cell = 0.5*(V/N)^1/3
Lchar2 = 0.5*power(abs(det(lat2))/sum(par_numIons2), 1/3);

Nslabs1 = round(L1/(Lchar1+(L1-Lchar1)*(cos(correlation_coefficient*pi/2))^2));
Nslabs2 = round(L2/(Lchar2+(L2-Lchar2)*(cos(correlation_coefficient*pi/2))^2));

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
    elseif cor_dir <= 0
        ord1 = 0;
    else
        ord1 = 1;
    end
    if length(whichIons2) > 0
        ord2 = sum(par_order2(whichIons2))/length(whichIons2);
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
            ord1_extra = sum(par_order1(whichIons1_extra))/length(whichIons1_extra);

        elseif cor_dir <= 0
            ord1_extra = 0;
        else
            ord1_extra = 1;
        end
        if ((ord1 > ord1_extra) && (cor_dir <= 0)) || ((ord1 < ord1_extra) && (cor_dir == 1))
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
            ord2_extra = sum(par_order2(whichIons2_extra))/length(whichIons2_extra);
        elseif cor_dir <= 0
            ord2_extra = 0;
        else
            ord2_extra = 1;
        end
        if ((ord2 > ord2_extra) && (cor_dir <= 0)) || ((ord2 < ord2_extra) && (cor_dir == 1))
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
 ionCount  = zeros(length(atomType),1);

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

 for ind = 1:length(atomType)
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
     ionCount(ind) = length(ions1) + length(ions2);
     numatom = ORG_STRUC.numIons(ind)*det(reshape(cell,[2,2]));
	if ionCount(ind) > numatom
          howmany = ionCount(ind) - numatom;
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
            ionCount(ind)=numatom;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fill the selected ions coordinates into the offspring
    addOn = cat(1,parent1(whichIons1(smaller1(ions1(:))) , :), parent2(whichIons2(smaller2(ions2(:))) , :));
    addOn = cat(1,addOn,candidates);
    offspring = cat(1,offspring,addOn);
     end
 end

numIons = transpose(ionCount);
potentialLattice=lat1;
end
