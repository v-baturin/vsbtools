function  dx_pro = NEB_trustRadius(dX, N)


global ORG_STRUC
global POP_STRUC


trustRadius = 0.25;     %% max displacement ~ 0.25A/atom
maxRelaxVolume = 10/100;   %%-dP/K  dP =2GPa K = 200~400GPa;

%--------     System Parameters      ------%

sumIons  = sum(ORG_STRUC.numIons);
dimension  = ORG_STRUC.dimension;
numDimension = 3*(sumIons+dimension);
dx_pro = zeros(numDimension,1);

%-------------------------------------------
dx_pro=dX;
tmp = ( reshape( dX',3,sumIons+dimension ) )';
CARTECOORDSDisplace = tmp(1+dimension:sumIons+dimension, :);
LATTICEStrain    = tmp(1:dimension,:);

%-------------------------------------------
% Cell Lattice Strain
% Relax_para = -dP/( K*trace(epsil) )
% dV/V ~ trace(epsil)  when dV is very small

epsil = LATTICEStrain;%/POP_STRUC.POPULATION(N).VOLUME;

if trace(epsil) > maxRelaxVolume
    dx_pro(1:9) = dX(1:9) * maxRelaxVolume/trace(epsil);
    %   disp('displacement of the lattic cell :')
    %   trace(epsil)
end


displacement = zeros(1, sumIons);
for i = 1:sumIons
    displacement(i) = norm( CARTECOORDSDisplace(i,:) ) ;
end

if max( displacement ) > trustRadius
    dx_pro(10:end) = dX(10:end) * trustRadius/max( displacement );
end

%-------------------------------------------


