function R = Operate_molecule(R0, Opt, lat)
%This routine is to apply the symmetry operations on the molcules R (in Cartesian Set)
%We need to do crystallographic operations always in fraction Set
%----No longer needs to consider the swap between a/b/c when it is monoclinic+orthorhombic.

%Step1: Obtain fraction axis
R0_frac = Cart2Frac(R0, lat);

%Step3: apply crystallographic operations 
for i = 1:size(R0,1)
   R_frac(i,:)  = R0_frac(i,:)*Opt(1:3,:) + Opt(4,:);
end

%Step5: obtain caretesian 
R   = Frac2Cart(R_frac, lat);
