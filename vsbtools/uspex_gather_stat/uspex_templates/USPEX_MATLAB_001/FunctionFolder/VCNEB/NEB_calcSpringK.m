function NEB_calcSpringK()

%--------------------------------------------%
global ORG_STRUC
global POP_STRUC

k_max=ORG_STRUC.K_max;
k_min=ORG_STRUC.K_min;

numImages = ORG_STRUC.numImages;
%--------------------------------------------%
E=zeros(numImages,1);
for i = 1:numImages
    E(i)= POP_STRUC.POPULATION(i).Enthalpy;
end
E_max = max(E);
E_min = min(E);


%-- Spring Constant Method  ---------------
%    1: Constant-K              (default)
%    2: Variable-K
%------------------------------------------

if numImages <= 2
    POP_STRUC.POPULATION(1).springK=0;
    POP_STRUC.POPULATION(numImages).springK=0;
else
    for i = 1:numImages
        E_n = POP_STRUC.POPULATION(i).Enthalpy;
        if 1==ORG_STRUC.optVarK
            POP_STRUC.POPULATION(i).springK = 0.5*( k_max+k_min - (k_max-k_min)*cos( pi*((E_n-E_min)/(E_max-E_min)) ) );
        else
            POP_STRUC.POPULATION(i).springK = 0.5*( k_max+k_min );
        end
    end
end
