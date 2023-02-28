function [coordinatonNumbers] = coordinationNumbersForComposition(numIons)

global USPEX_STRUC
global ORG_STRUC

if isempty(ORG_STRUC.maxAt)
    grid_sz = sum(ORG_STRUC.numIons)+1;
else
    grid_sz = ORG_STRUC.maxAt+1;
end

switch size(numIons,2)
    case 5
        cn_array = reshape([USPEX_STRUC.EXT_RandTop.CN_grid{:,:}],grid_sz,grid_sz,grid_sz,grid_sz,grid_sz,grid_sz,[]);
        coordinatonNumbers = cn_array(numIons(1)+1,numIons(2)+1,numIons(3)+1,numIons(4)+1,numIons(5)+1,:);
    case 4
        cn_array = reshape([USPEX_STRUC.EXT_RandTop.CN_grid{:,:}],grid_sz,grid_sz,grid_sz,grid_sz,[]);
        coordinatonNumbers = cn_array(numIons(1)+1,numIons(2)+1,numIons(3)+1,numIons(4)+1,:);
    case 3
        cn_array = reshape([USPEX_STRUC.EXT_RandTop.CN_grid{:,:}],grid_sz,grid_sz,grid_sz,[]);
        coordinatonNumbers = cn_array(numIons(1)+1,numIons(2)+1,numIons(3)+1,:);
    case 2
        cn_array = reshape([USPEX_STRUC.EXT_RandTop.CN_grid{:,:}],grid_sz,grid_sz,[]);
        coordinatonNumbers = cn_array(numIons(1)+1,numIons(2)+1,:);
    case 1
        cn_array = reshape([USPEX_STRUC.EXT_RandTop.CN_grid{:,:}],grid_sz,[]);
        coordinatonNumbers = cn_array(numIons(1)+1,:);
end
