function [coor, numIons] = RemoveDup(coor, numIons, tolerance)
bad_rank = [];
for i = 2:size(coor,1)
    for j = 1:i-1
        dist = coor(i,:)-coor(j,:);
        dist = dist - round(dist);
        if (norm(dist)<tolerance)
           bad_rank = [bad_rank;i];
           ID = findatomType(i, numIons);
           numIons(ID) = numIons(ID)-1;
           break;
        end
    end
end
coor(bad_rank,:)=[];

