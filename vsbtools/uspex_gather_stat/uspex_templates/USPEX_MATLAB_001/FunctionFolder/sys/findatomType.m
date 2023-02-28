function type = findatomType(id, numIons)

type = 1;
SUM = 0;
for i = 1:length(numIons)
    numIons(i)=numIons(i)+SUM;
    if (id > SUM) && (id <=numIons(i))
        type = i;
        break
    end
    SUM = numIons(i);
end


