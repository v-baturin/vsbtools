function goodComposition = checkComposition(numIons)

global ORG_STRUC

goodComposition = 1;
for i = 1 : length(ORG_STRUC.atomType)
    if (numIons(i) < ORG_STRUC.numIons(1,i)) || (numIons(i) > ORG_STRUC.numIons(2,i))
        goodComposition = 0;
        break;
    end
end