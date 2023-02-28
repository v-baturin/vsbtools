function coef = find_coef_antiseeds(numIons, fitness, ranking)

%finding coefficients of antiSeedsMax for every composition

global ORG_STRUC
global POP_STRUC

num = 0;
for i = 1 : round(length(fitness)*ORG_STRUC.bestFrac)
    if isequal(numIons, POP_STRUC.POPULATION(ranking(i)).numIons)
        num = num + 1;
    end
end

%finding number of all possible compositions
numIons_range = ORG_STRUC.numIons;

if length(ORG_STRUC.atomType) == 1
    num_comp = numIons_range(2,1)-numIons_range(1,1)+1;
elseif length(ORG_STRUC.atomType) == 2
    num_comp = (numIons_range(2,1)-numIons_range(1,1)+1)*(numIons_range(2,2)-numIons_range(1,2)+1);
elseif length(ORG_STRUC.atomType) == 3
    num_comp = (numIons_range(2,1)-numIons_range(1,1)+1)*(numIons_range(2,2)-numIons_range(1,2)+1)*(numIons_range(2,3)-numIons_range(1,3)+1);
end

num1 = round(length(fitness)*ORG_STRUC.bestFrac)/num_comp; % average structures per one composition
coef = num/num1;
