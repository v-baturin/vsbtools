function [whichStructure, whichAtom] = find_candidate_remove()

global AR_VARIANTS

%finding the number of structure to choose as parent
delta = [];
for i = 1 : length(AR_VARIANTS.Structure)
    delta(i)=0;
    for j = 1 : length(AR_VARIANTS.Structure(i).Atom)
        delta(i) = delta(i) + AR_VARIANTS.Structure(i).Atom(j).remove;
    end
end    

if sum(delta)==0 %no more variants for removing atom
    whichStructure = -1;
    whichAtom = -1;
else
    rand_numb = rand(1)*sum(delta);
    
    for i = 1 : length(AR_VARIANTS.Structure)
        if rand_numb <= sum(delta(1:i))
            whichStructure = i;
            break;
        end
    end
    
    %finding the number of atom to delete
    rand_numb1 = rand_numb - sum(delta(1:i-1));
    
    delta1 = [];
    for j = 1 : length(AR_VARIANTS.Structure(i).Atom)
        delta1(j) = AR_VARIANTS.Structure(i).Atom(j).remove;
        if rand_numb1 <= sum(delta1(1:j))
            whichAtom = j;
            break;
        end
    end
end