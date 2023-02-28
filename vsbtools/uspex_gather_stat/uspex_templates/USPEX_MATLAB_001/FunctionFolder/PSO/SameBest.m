function same_in=SameBest(fitness)
%This program is to find the best structures has not changed for how many generations
%Last updated by Qiang Zhu (2013/10/24)
global POP_STRUC
global ORG_STRUC

same_in = 0;
[tmp,ind] = min(fitness);
f1 = POP_STRUC.POPULATION(ind).FINGERPRINT;
for i = POP_STRUC.generation - ORG_STRUC.stopCrit + 1 : POP_STRUC.generation - 1
    try
        load([POP_STRUC.resFolder '/generation' num2str(i) '/POP_STRUC.mat']);
        f2 = POP_STRUC.POPULATION(POP_STRUC.ranking(1)).FINGERPRINT;
        if isempty(f1) | isempty(f2)
            cos_dist = 1;
        else
            cos_dist = cosineDistance(f1, f2, ORG_STRUC.weight);
        end
        if (cos_dist < ORG_STRUC.toleranceFing)
            same_in = same_in + 1;
        end
    catch
    end
end
load Current_POP.mat

