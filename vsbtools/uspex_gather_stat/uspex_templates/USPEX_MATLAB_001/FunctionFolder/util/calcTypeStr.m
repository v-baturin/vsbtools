function calcType_str = calcTypeStr()

global ORG_STRUC

calcType_str = [num2str(ORG_STRUC.dimension) '' num2str(ORG_STRUC.molecule) '' num2str(ORG_STRUC.varcomp)];
if str2num(calcType_str) < 0
    calcType = -1*str2num(calcType_str);
    calcType_str = ['M' num2str(calcType)];
end
