function [whichStructure, whichAtom, whichAtom2, whichEdge, typeAtom_add] = find_candidate_add()

global ORG_STRUC
global AR_VARIANTS

%finding the number of structure to choose as parent
delta = [];
for i = 1 : length(AR_VARIANTS.Structure)
    delta(i)=0;
    for j = 1 : length(AR_VARIANTS.Structure(i).Atom)
        for k = 1 : length(AR_VARIANTS.Structure(i).Atom(j).Edge)
            delta(i) = delta(i) + sum(AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(:));
        end
    end
end

if sum(delta)==0 %no more variants for adding atom
    whichStructure = -1;
    whichAtom = -1;
    whichAtom2 = -1;
    whichEdge = -1;
    typeAtom_add = -1;
else
    rand_numb = rand(1)*sum(delta);
    
    for i = 1 : length(AR_VARIANTS.Structure)
        if rand_numb <= sum(delta(1:i))
            whichStructure = i;
            break;
        end
    end
    
    %finding the number of atom to add new atom near it
    rand_numb1 = rand_numb - sum(delta(1:i-1));
    
    delta1 = [];
    for j = 1 : length(AR_VARIANTS.Structure(i).Atom)
        delta1(j) = 0;
        for k = 1 : length(AR_VARIANTS.Structure(i).Atom(j).Edge)
            delta1(j) = delta1(j) + sum(AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(:));
        end
        if rand_numb1 <= sum(delta1(1:j))
            whichAtom = j;
            break;
        end
    end
    
    %finding the edge to add new atom near its center
    rand_numb2 = rand_numb1 - sum(delta1(1:j-1));
    
    delta2 = [];
    for k = 1 : length(AR_VARIANTS.Structure(i).Atom(j).Edge)
        delta2(k) = sum(AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(:));
        if rand_numb2 <= sum(delta2(1:k))
            whichAtom2 = AR_VARIANTS.Structure(i).Atom(j).Edge(k).Atom;
            whichEdge = k;
            break;
        end
    end
    
    %finding the type of atom to add
    rand_numb3 = rand_numb2 - sum(delta2(1:k-1));
    
    delta3 = [];
    for l = 1 : length(ORG_STRUC.atomType)
        delta3(l) = AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(l);
        if rand_numb3 <= sum(delta3(1:l))
            typeAtom_add = l;
            break;
        end
    end
end