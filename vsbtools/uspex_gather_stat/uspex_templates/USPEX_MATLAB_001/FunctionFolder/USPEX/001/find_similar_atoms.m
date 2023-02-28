function sim_atoms = find_similar_atoms(Structure, whichAtom)

global AR_VARIANTS

numIons = AR_VARIANTS.Structure(Structure).numIons;
coord   = AR_VARIANTS.Structure(Structure).coord;
lattice = AR_VARIANTS.Structure(Structure).lattice;
typeAtom = find_atomType(numIons, whichAtom);

sim_atoms = [];
ns = 0;

for curAtom = sum(numIons(1:typeAtom-1))+1 : sum(numIons(1:typeAtom))
    d1 = [];
    d2 = [];
    for type = 1 : length(numIons)
        n = 0;
        for i = sum(numIons(1:type-1))+1 : sum(numIons(1:type))
            n = n + 1;
            d1(n) = distAtoms(coord(whichAtom,:)*lattice, coord(i,:)*lattice);
            d2(n) = distAtoms(coord(curAtom,:)*lattice, coord(i,:)*lattice);            
        end
        diff = sort(d1)-sort(d2);
        dist_fp(type) = norm(diff)/length(diff);
    end
    
    if sum(dist_fp) < 0.01 %similar atoms
        ns = ns + 1;
        sim_atoms(ns) = curAtom;
    end
end