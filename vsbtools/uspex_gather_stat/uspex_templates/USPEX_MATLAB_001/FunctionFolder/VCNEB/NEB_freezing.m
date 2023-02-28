function NEB_freezing(step)

%--------------------------------------------------------
global ORG_STRUC
global POP_STRUC
%--------     System Parameters      --------%
numImages = ORG_STRUC.numImages;
%--------------------------------------------%
for i=1:numImages
    convThreshold = max([POP_STRUC.POPULATION(i).errcFnebMax, POP_STRUC.POPULATION(i).erraFnebMax,...
        POP_STRUC.POPULATION(i).errcFnebRms, POP_STRUC.POPULATION(i).erraFnebRms]);
    if (convThreshold < ORG_STRUC.ConvThreshold)
        if (ORG_STRUC.optFreezing==1) || (ORG_STRUC.CalcType==2)
            POP_STRUC.POPULATION(i).freezing = 1;
        else
            POP_STRUC.POPULATION(i).freezing = 0;
        end
    else
        POP_STRUC.POPULATION(i).freezing = 0;
    end
end
return;
