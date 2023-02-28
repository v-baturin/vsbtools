function CreateCalcFolder(varargin)

global ORG_STRUC


cd (ORG_STRUC.homePath)


if length(varargin) > 0
    mkdir(['CalcFoldTemp']);
    mkdir(['CalcFoldTemp/.data']);
    unixCmd(['cp ' ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup/data/*  CalcFoldTemp/.data/']);
else
    for folds = 1 : ORG_STRUC.numParallelCalcs
        if ~exist(['CalcFold' num2str(folds)])
			%%%%%%%%%%%%%%
            mkdir(['CalcFold' num2str(folds)]);
            unixCmd(['cp  ', ORG_STRUC.USPEXPath, '/FunctionFolder/Tool/getStuff CalcFold' num2str(folds)]);
			unixCmd(['chmod +x CalcFold' num2str(folds) '/getStuff ']);
			%%%%%%%%%%%%%%
            if sum(ORG_STRUC.abinitioCode==1) == 0 || size(ORG_STRUC.abinitioCode,2) ~= sum(ORG_STRUC.abinitioCode==1)  % case of mixed calculators
                unixCmd(['cp ./'  ORG_STRUC.specificFolder '/* CalcFold' num2str(folds)]);
            else
                %% for VASP code, only copy the vdw_kernel.bindat file
                if exist([ORG_STRUC.specificFolder '/vdw_kernel.bindat'])
                    unixCmd(['cp ./'  ORG_STRUC.specificFolder '/vdw_kernel.bindat CalcFold' num2str(folds)]);
                end
            end
			%%%%%%%%%%%%%%
        end
    end
end



