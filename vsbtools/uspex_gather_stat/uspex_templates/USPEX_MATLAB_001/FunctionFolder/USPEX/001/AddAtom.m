function [goodComposition, son, numIons_son, LATTICE] = AddAtom(whichStructure, whichAtom, whichAtom2, typeAtom_add)

global ORG_STRUC
global AR_VARIANTS

father         = AR_VARIANTS.Structure(whichStructure).coord;
numIons_father = AR_VARIANTS.Structure(whichStructure).numIons;
LATTICE        = AR_VARIANTS.Structure(whichStructure).lattice;

typeAtom = find_atomType(numIons_father, whichAtom);
numIons_son = numIons_father;
numIons_son(typeAtom_add) = numIons_father(typeAtom_add) + 1;

son = [];
goodComposition = checkComposition(numIons_son);

if goodComposition == 1
    covRad = [];
    for i = 1 : length(ORG_STRUC.atomType)
        covRad(i) = str2num(covalentRadius(ORG_STRUC.atomType(i)));
    end
    
    %we add atom to the center of the edge between whichAtom and whichAtom2 atoms
    typeAtom  = find_atomType(numIons_father,whichAtom);
    typeAtom2 = find_atomType(numIons_father,whichAtom2);
        
    father_cart = father*LATTICE;
    for i = 1:3
        center_mass(i) = sum(father_cart(:,i))/sum(numIons_father);
    end
    center_edge = (father_cart(whichAtom2,:) + father_cart(whichAtom,:))/2;
    ddd = center_edge - center_mass;
    
    edge_vector1 = father_cart(whichAtom2,:) - father_cart(whichAtom,:);
    edge_vector = edge_vector1/norm(edge_vector1);
    
    if norm(cross(ddd,edge_vector)) == 0
        ddd = [rand(1), rand(1), rand(1)];
    end
    
    nnn = ddd - (ddd*edge_vector')*edge_vector;
    nnn = nnn/norm(nnn);
    
    rrr = covRad(typeAtom_add) + max(covRad(typeAtom), covRad(typeAtom2));
    if (rrr > norm(edge_vector1)/2)
        lll = sqrt(rrr^2 - (norm(edge_vector1)/2)^2);
    else
        lll = 0;
    end
        
    adding_atom = (center_edge + lll*nnn + 0.01*[rand(1),rand(1),rand(1)]) * LATTICE^(-1);
    %adding_atom = (center_edge + lll*nnn) * LATTICE^(-1);
    
    % find coordinates of structure with additional atom
    for i = 1 : sum(numIons_father(1:typeAtom_add))
        son(i,:) = father(i,:);
    end
    son(end+1,:) = adding_atom;
    
    for i = sum(numIons_father(1:typeAtom_add))+1 : sum(numIons_father)
        son(i+1,:) = father(i,:);
    end
    
    if ~isreal(son)
        disp('complex coordinates!!!!');
    end
end