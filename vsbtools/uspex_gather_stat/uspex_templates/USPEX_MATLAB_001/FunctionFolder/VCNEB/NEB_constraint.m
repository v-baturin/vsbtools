function NEB_constraint()

%-------------------------------------%
global ORG_STRUC

sumIons = sum(ORG_STRUC.numIons);
dimension = ORG_STRUC.dimension;

numDimension=3*(sumIons+dimension);

%-------------------------------------%
optRelaxType=ORG_STRUC.optRelaxType;
%-------------------------------------%
ORG_STRUC.whetherConstraint = zeros(numDimension,1);
%--------------------------------------------------------------------
%Force = [ Cell Force *9, Natom Force * 3sumIons]
%---[  1  2  3  ]
%---|  4  5  6  ]
%---[  7  8  9  ]

if     1==optRelaxType         % atom relax
    ORG_STRUC.whetherConstraint( 3*dimension+1:3*(sumIons+dimension) ) = 1;
elseif 2==optRelaxType         % cell relax
    ORG_STRUC.whetherConstraint( 1:3*dimension ) = 1;
elseif 3==optRelaxType         % atom+cell relax
    ORG_STRUC.whetherConstraint( : )=1;
end

%---------------------------------------------------------------------
