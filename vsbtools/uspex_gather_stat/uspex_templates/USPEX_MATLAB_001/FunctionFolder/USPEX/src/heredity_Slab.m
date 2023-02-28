function [parent, whichIons] = heredity_Slab(parent, order, Nslab, fracFrac, dimension, cor_dir, ID)
% chooses the ordered slab for both parents from (anti)correlation

% The order from intail shift
if ID==1
   whichIons = find(parent(:,dimension)<fracFrac);
else
   whichIons = find(parent(:,dimension)>fracFrac);
end
ord = heredity_CalcOrder(whichIons, order, cor_dir);

% Here we try different slabs, if the resulting order is better, use this slab
for i = 1 : Nslab
    offset_extra = rand;
    parent(:,dimension) = parent(:,dimension) + offset_extra;
    parent(:,dimension) = parent(:,dimension) - floor(parent(:,dimension));
    if ID==1
       whichIons_extra = find(parent(:,dimension)<fracFrac);
    else
       whichIons_extra = find(parent(:,dimension)>fracFrac);
    end
   
    ord_extra = heredity_CalcOrder(whichIons_extra, order, cor_dir);

    if ((ord > ord_extra) && (cor_dir <= 0)) || ((ord < ord_extra) && (cor_dir == 1))
        parent(:,dimension) = parent(:,dimension) - offset_extra;
        parent(:,dimension) = parent(:,dimension) - floor(parent(:,dimension));
    else
        whichIons = whichIons_extra;
        ord1 = ord_extra;
    end
end

%-----------------------------------------------------------
function ord = heredity_CalcOrder(whichIons, order, cor_dir)
%-----------------------------------------------------------
%small function to quickly calculate order

if length(whichIons) > 0
    ord = sum(order(whichIons))/length(whichIons);
elseif cor_dir <= 0
    ord = 0;
else
    ord = 1;
end

