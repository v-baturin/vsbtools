function Coords= NEW_ZMATRIXCOORD(ZMATRIX,format)

%This function is to transfor Zmatrix to XYZ
%Remember that Zmatrix is a define in spherical coordinates systems
%So we need a lot of transformation from (r,theta,phi) to (x,y,z)

ZMATRIX = real(ZMATRIX);
N_atom  = size(ZMATRIX, 1);
Coords  = zeros(N_atom, 3);
origin  = ZMATRIX(1,:);
if N_atom > 1
   Coords(2,3) = ZMATRIX(2,1)*cos(ZMATRIX(2,2));
   Coords(2,1) = ZMATRIX(2,1)*sin(ZMATRIX(2,2))*cos(ZMATRIX(2,3));
   Coords(2,2) = ZMATRIX(2,1)*sin(ZMATRIX(2,2))*sin(ZMATRIX(2,3));
   if N_atom > 2
      for i = 3:size(ZMATRIX,1)
          if i==3
             Ref = Coords(format(3,1:2),:);
          else
             Ref = Coords(format(i,:),:);
          end
          Coords(i, :)= GetXYZ(Ref, ZMATRIX(i,:));
      end
   end   
end
Coords = bsxfun(@plus, real(Coords), origin);

%--------------------------------------------------------
function coor = GetXYZ(Ref,Zmatrix)
%--------------------------------------------------------
%from reference (i = 1, 2, 3), we construct the followings
%1, origin (1)
%2, z-axis along 1->2
%3, y-axis is perpendicular to plane (1-2-3)
%the new sperical system as (r, theta, phi)
r     = Zmatrix(1);
theta = Zmatrix(2);
phi   = -1*Zmatrix(3);

coor = [r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)];
U1 = Ref(2,:) - Ref(1,:);
if size(Ref,1)==2
   U2 = [1 0 0];
else
   U2 = Ref(3,:) - Ref(2,:);
end
Z  = U1/norm(U1);
Y  = cross(U1,U2);
Y  = Y/norm(Y);
X  = cross(Y,Z);
X  = X/norm(X);
coor = coor/([X;Y;Z]');
coor = coor + Ref(1,:);


