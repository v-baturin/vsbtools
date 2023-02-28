function weight = CalcWeight_001(numIons)

L = length(numIons);
S = 0;
weight = zeros(L*L,1);
for i = 1:L
    for j = 1:L
        ind = (i-1)*L+j;
        weight(ind) = (numIons(i)*numIons(j));
        S = S + (numIons(i)*numIons(j));
    end
end
    
weight = weight/S;