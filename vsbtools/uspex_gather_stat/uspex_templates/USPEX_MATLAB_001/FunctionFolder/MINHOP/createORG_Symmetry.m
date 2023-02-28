function createORG_Symmetry(inputFile)

% USPEX Version 9.3.0
% new tags: dimension/varcomp/molecule
% deleted tags: supercomputer name
global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%[nothing, doSpaceGroup] = unix (['./getStuff ' inputFile ' doSpaceGroup 1']);
doSpaceGroup = python_uspex(getPy, ['-f ' inputFile ' -b doSpaceGroup -c 1']);
if isempty(doSpaceGroup)
    if (ORG_STRUC.dimension==0) || (ORG_STRUC.dimension==2) || (ORG_STRUC.dimension==-3)
        doSpaceGroup = '0'; % default
    else
        doSpaceGroup = '1'; % default
    end
end
ORG_STRUC.doSpaceGroup = str2num(doSpaceGroup);


% do we want to symmetrize the structure using Stokes SG determination code (using symmetrized.cif)
% done at the last optimization step (when not shaking - symmetrize)
%[nothing, symmetrize] = unix(['./getStuff ' inputFile ' symmetrize 1']);
symmetrize = python_uspex(getPy, ['-f ' inputFile ' -b symmetrize -c 1']);
if ~isempty(symmetrize)
    ORG_STRUC.symmetrize = str2num(symmetrize);
end
% Space group determination tolerance
%[noathing, SGtolerance] = unix(['./getStuff ' inputFile ' SymTolerance 1']);
SGtolerance = python_uspex(getPy, ['-f ' inputFile ' -b SymTolerance -c 1']);
if ~isempty(SGtolerance)
    SGtolerance = deblank(SGtolerance);
    if strcmp(lower(SGtolerance), 'high')
        ORG_STRUC.SGtolerance = 0.05;
    elseif strcmp(lower(SGtolerance), 'medium')
        ORG_STRUC.SGtolerance = 0.10;
    elseif strcmp(lower(SGtolerance), 'low')
        ORG_STRUC.SGtolerance = 0.20;
    else % number!
        ORG_STRUC.SGtolerance = str2num(SGtolerance);
    end
end
if isempty(ORG_STRUC.SGtolerance) %In case some illegal strings
    ORG_STRUC.SGtolerance = 0.10;
end
% coefficient between mindist and symmetrization distance (1 by default, sometimes > 1 needed)
%[nothing, sym_coef] = unix(['./getStuff ' inputFile ' constraint_enhancement 1']);
sym_coef = python_uspex(getPy, ['-f ' inputFile ' -b constraint_enhancement -c 1']);
if ~isempty(sym_coef)
    ORG_STRUC.sym_coef = str2num(sym_coef);
end

