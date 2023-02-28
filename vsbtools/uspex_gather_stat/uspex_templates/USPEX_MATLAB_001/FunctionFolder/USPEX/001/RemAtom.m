function [goodComposition, son, numIons_son, LATTICE] = RemAtom(whichStructure, whichAtom)

global AR_VARIANTS

father         = AR_VARIANTS.Structure(whichStructure).coord;
numIons_father = AR_VARIANTS.Structure(whichStructure).numIons;
LATTICE        = AR_VARIANTS.Structure(whichStructure).lattice;

typeAtom = find_atomType(numIons_father, whichAtom);
numIons_son = numIons_father;
numIons_son(typeAtom) = numIons_father(typeAtom) - 1;

son = [];
goodComposition = checkComposition(numIons_son);
    
if goodComposition == 1
    for i = 1 : whichAtom-1
        son(i,:) = father(i,:);
    end
    
    for i = whichAtom+1 : sum(numIons_father)
        son(i-1,:) = father(i,:);
    end
end