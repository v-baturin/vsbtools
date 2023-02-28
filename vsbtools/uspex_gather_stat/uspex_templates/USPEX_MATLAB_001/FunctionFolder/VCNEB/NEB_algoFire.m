function   dx_QMIN = NEB_algoFire(N)

%-----------------------------------------%
global ORG_STRUC
global POP_STRUC
%-----------------------------------------%

Nmin_FIRE = 10 ;     %- Fire Parameter :  Nmin=5;
finc_FIRE = 1.1;
fdec_FIRE = 0.5;
a0_FIRE = 0.1;
fa_FIRE = 0.99;

%--------     System Parameters      ------%
numImages = ORG_STRUC.numImages;
numIons  = ORG_STRUC.numIons;
sumIons  = sum(ORG_STRUC.numIons);
dimension = ORG_STRUC.dimension;
numDimension = 3*( sum(numIons) + dimension );

%---------------------------------------


%disp('FIRE -- : ');
%POP_STRUC.POPULATION(N).Fire

E=zeros(numImages,1);
for i = 1:numImages
    E(i)= POP_STRUC.POPULATION(i).Enthalpy;
end

%E_max = max(E);  E_min = min(E);
F = POP_STRUC.POPULATION(N).F_pro + POP_STRUC.POPULATION(N).F_ela;


if abs(POP_STRUC.POPULATION(N).Fire.dt)<0.001  % Prevent dt=0, no relaxation
    POP_STRUC.POPULATION(N).Fire.dt=1e-2;
end

%---------------------------------------------------------------------
if 0==POP_STRUC.step
    Ncount=0;
    a = a0_FIRE;
    P = 0;
    dt = ORG_STRUC.dt;
    v = zeros(numDimension,1);
else
    Ncount=POP_STRUC.POPULATION(N).Fire.Ncount;
    dt=POP_STRUC.POPULATION(N).Fire.dt;
    a=POP_STRUC.POPULATION(N).Fire.a;
    v=POP_STRUC.POPULATION(N).Fire.v;
    
    P = F'*v;
    v= (1-a)*v + a*F/norm(F)*norm(v);
    
    if P > 0
        Ncount= Ncount+1;
        if Ncount > Nmin_FIRE
            dt = min(dt*finc_FIRE, ORG_STRUC.dt*5);
            a = a*fa_FIRE;
        end
    else
        Ncount = 0;
        v(:)  = 0;
        dt = dt*fdec_FIRE;
        a = a0_FIRE;
    end
end

v = v + F*dt/2;
dx_QMIN = v*dt + F*dt*dt/2;


%%-------
POP_STRUC.POPULATION(N).Fire.Ncount=Ncount;
POP_STRUC.POPULATION(N).Fire.dt = dt;
POP_STRUC.POPULATION(N).Fire.a = a;
POP_STRUC.POPULATION(N).Fire.v = v;
POP_STRUC.POPULATION(N).Fire.P = P;

%POP_STRUC.POPULATION(N).Fire
