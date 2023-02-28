function Zmatrix = NEW_coord2Zmatrix(coords, format)

coords  = real(coords);

Zmatrix = coords;  %1st atom always = coords
N_atom  = size(coords,1);

if N_atom > 1
   coords = bsxfun(@minus, coords, coords(1,:));
%2nd atom, define it in spherical coordinates systems
   Zmatrix(2,1) = real(norm(coords(2,:)));
   if coords(2,3)==0
       Zmatrix(2,2)=pi/2;
   else
       Zmatrix(2,2)=acos(coords(2,3)/Zmatrix(2,1));
   end

   if coords(2,2)==0
       Zmatrix(2,3)=0;
   else
       Zmatrix(2,3)=atan2(coords(2,2),coords(2,1));
   end

   for ind = 3:N_atom
       a1 = coords(ind,:);
       a2 = coords(format(ind,1),:);
       a3 = coords(format(ind,2),:);
       Zmatrix(ind,1) = real(norm(a2-a1));
       Zmatrix(ind,2) = GetAngle(a1,a2,a3);
       if ind == 3 % the dihedral angle between 1-2-3 and XY plane
          a4 = a3 + [ 1 0 0];
          %Zmatrix(ind,3) = -1*GetDihedral(a1,a2,a3,a4);
       else
          a4 = coords(format(ind,3),:);
          %Zmatrix(ind,3) = -1*GetDihedral(a1,a2,a3,a4);
       end
          Zmatrix(ind,3) = GetDihedral(a1,a2,a3,a4);
   end
end
Zmatrix = real(Zmatrix);

%----------Bond Angle
function angle = GetAngle(a1,a2,a3)
v1=a1-a2;
v2=a3-a2;
angle = acos(dot(v1,v2)/norm(v1)/norm(v2));

%-----------Torsion
function torsion = GetDihedral(a1,a2,a3,a4)
p = a2-a1;
q = a3-a2;
r = a4-a3;
n1 = cross(p,q);
n2 = cross(q,r);
torsion = acos(dot(n1,n2)/(norm(n1)*norm(n2)));
center =(a1+a2+a3)/3;
if dot(n1, a4-center)<0
   torsion = torsion*-1;
end

