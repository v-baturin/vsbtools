function typeAtom = find_atomType(numIons, whichAtom)

for i = 1 : length(numIons)
    if whichAtom <= sum(numIons(1:i))
        typeAtom = i;
        break;
    end
end