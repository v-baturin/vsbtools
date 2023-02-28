function [goodSon, transNow, numIons_son, LATTICE] = swapIons_transmutation_final_001(Ind_No)

% USPEX Version 7.2.4
% ORG_STRUC.maxIons and ORG_STRUC.minIons are not activated yet, not sure if they are needed at all
% the algorithm: 
% 1. Determine how many blocks are removed and how many blocks are added
% 2. Determine which atoms are being transmutated and make them 'grey' (atom Type = 0)
% 3. Assign types to grey atoms, if there is excess of them - delete, if we need extra atoms to add - add them randomly
global POP_STRUC
global ORG_STRUC

numIons = ORG_STRUC.numIons;
%addBlocks = zeros(1, length(numBlocks));  % how many blocks to add
%delBlocks = addBlocks;                    % how many blocks to remove
father = POP_STRUC.POPULATION(Ind_No).COORDINATES;
%goodFather = father;
numIons_father = POP_STRUC.POPULATION(Ind_No).numIons;
LATTICE = POP_STRUC.POPULATION(Ind_No).LATTICE;

% we set the max fraction of atoms to be transmutated, not the set number of it, thus we calculate the set number first
% howManyTrans == 0 means the user doesn't want to think for himself (we told him so in INPUT.txt) 

%howManyTrans = round(ORG_STRUC.howManyTrans*sum(numIons0));
%if howManyTrans < 1 | isnan( howManyTrans )
%    howManyTrans = 1;
%end

howManyAtomsTrans = round(0.2*sum(numIons_father));

% would be kind of boring always to swap the same amount of ions, so we introduce an element of randomness
transNow = round(rand(1)*(howManyAtomsTrans-1)+1);

% sometimes transition doesn't exist! f.e.: 5 0 0 and only allowed transmutation is 2 <=> 3
% check whether structure can be transmuted, i.e. exist non zero block that is listed in allowed transmutations

if length(ORG_STRUC.atomType) == 2 %old code
    
    count = 1;
    goodComp = 0;
    while (count < 10000) && (goodComp == 0)
        count = count + 1;
        atoms2mutate = randperm(sum(numIons_father));
        num_from1_to2 = 0;
        num_from2_to1 = 0;
        atoms_from1_to2 = [];
        atoms_from2_to1 = [];
        for i = 1:transNow
            if atoms2mutate(i) <= numIons_father(1)
                num_from1_to2 = num_from1_to2 + 1;
                atoms_from1_to2(num_from1_to2) = atoms2mutate(i);
            else
                num_from2_to1 = num_from2_to1 + 1;
                atoms_from2_to1(num_from2_to1) = atoms2mutate(i);
            end
        end
        numIons_son(1) = numIons_father(1) - num_from1_to2 + num_from2_to1;
        numIons_son(2) = numIons_father(2) + num_from1_to2 - num_from2_to1;
        
        if size(ORG_STRUC.numIons,2) == 2
            if (numIons_son(1)>=numIons(1,1)) && (numIons_son(1)<=numIons(2,1)) && (numIons_son(2)>=numIons(1,2)) && (numIons_son(2)<=numIons(2,2))
                goodComp = 1;
            end
        elseif size(ORG_STRUC.numIons,2) == 3
            if (numIons_son(1)>=numIons(1,1)) && (numIons_son(1)<=numIons(2,1)) && ...
                    (numIons_son(2)>=numIons(1,2)) && (numIons_son(2)<=numIons(2,2)) && ...
                    (numIons_son(3)>=numIons(1,3)) && (numIons_son(3)<=numIons(2,3))
                goodComp = 1;
            end
        end
    end
    %disp('count');
    %disp(count);
    
    %disp('transNow');
    %disp(transNow);
    
    %disp('n_type1_father');
    %disp(numIons_father(1));
    
    %disp('n_type2_son');
    %disp(numIons_father(2));
    
    %disp('atoms from 1 to 2; from 2 to 1');
    %disp(atoms_from1_to2);
    %disp(atoms_from2_to1);
    
    %disp('n_type1_son');
    %disp(numIons_son(1));
    
    %disp('n_type2_son');
    %disp(numIons_son(2));
    
    goodSon = [];
    transNow = 0;
    if (goodComp == 1) && (count<10000)
        %numions_son(1)
        j = 0;
        for i = 1:numIons_father(1)
            if isempty(find(atoms_from1_to2==i))
                j = j + 1;
                goodSon(j,:) = father(i,:);
            end
        end
        for i = 1:num_from2_to1
            j = j + 1;
            goodSon(j,:) = father(atoms_from2_to1(i),:);
        end
        
        %disp('n_type1__');
        %disp(j);
        %numions_son(2)
        for i = 1:numIons_father(2)
            if isempty(find(atoms_from2_to1==numIons_father(1)+i))
                j = j + 1;
                goodSon(j,:) = father(numIons_father(1)+i,:);
            end
        end
        for i = 1:num_from1_to2
            j = j + 1;
            goodSon(j,:) = father(atoms_from1_to2(i),:);
        end
        %disp('n_type_sum__');
        %disp(j);
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(ORG_STRUC.atomType) >= 3 %new code - universal
    possible_trans = [0, 1, 0; ...
                      1, 0, 0; ...
                      1, 1, 0];
    type_atom = [];
    for i = 1:sum(numIons_father)
        type_atom(i) = find_atomType(numIons_father, i);
    end
    
    for i = 1:1000 %number attempts to choose atoms to transmutate
        %generate set of atoms to be transmutated
        perm = randperm(sum(numIons_father));
        atoms_trans = [];
        for ii = 1:transNow
            atoms_trans(ii) = perm(ii);
        end
        
        %2. try all possible transmutations for the given set of atoms
        for j = 1:100
            % choose random transmutation for the given set of atoms
            how_trans = [];
            for ii = 1:transNow
                how_trans(ii) = 1+round(rand(1)*(length(ORG_STRUC.atomType)-1));
            end
            
            goodComp = 1;
            numIons_son = numIons_father;
            type_atom_new = type_atom;
            
            %2.1 check if this transmutation set is possible
            for k = 1:transNow
                type_initial = type_atom(atoms_trans(k)); %atom type before transmutation
                type_final = how_trans(k); %atom type after transmutation
                if (type_initial == type_final) || ~possible_trans(type_initial,type_final)
                    goodComp = 0;
                    break;
                else
                    numIons_son(type_initial) = numIons_son(type_initial) - 1;
                    numIons_son(type_final) = numIons_son(type_final) + 1;
                    type_atom_new(atoms_trans(k)) = type_final;
                end
            end
            
            if goodComp
                goodComp = checkComposition(numIons_son);
            end
            
            if goodComp
                break
            end            
        end
        
        if goodComp
            break
        end
    end
    
    goodSon = [];
    if goodComp
        %using new type_atom we create goodSon
        disp(['Ind_No: ' num2str(Ind_No)]);
        disp(['numIons_father: ' num2str(numIons_father)]);
        disp(['transNow: ' num2str(transNow)]);
        disp(['type_atom: ' num2str(type_atom)]);
        disp(['type_atom_new: ' num2str(type_atom_new)]);
        [nothing, new_order] =  sort(type_atom_new);
        disp(['new_order: ' num2str(new_order)]);
        goodSon = [];
        for i = 1:length(ORG_STRUC.atomType)
            numIons_son_(i) = 0;
        end
        for i = 1:sum(numIons_father)
            goodSon(i,:) = father(new_order(i),:);
            numIons_son_(type_atom_new(new_order(i))) = numIons_son_(type_atom_new(new_order(i))) + 1;
        end
        disp(['numIons_son: ' num2str(numIons_son)]);
        disp(['numIons_son_: ' num2str(numIons_son_)]);
    end
end
