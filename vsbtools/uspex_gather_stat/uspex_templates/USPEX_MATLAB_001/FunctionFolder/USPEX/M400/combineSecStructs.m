function sec_parts = combineSecStructs(res_num)
% $Rev$
% $Author$
% $Date$

while 1
    a = RandInt(1,7,[0 res_num]);
    sec_parts = round(a/sum(a)*res_num);
    if sum(sec_parts) == res_num
        break
    end
end
