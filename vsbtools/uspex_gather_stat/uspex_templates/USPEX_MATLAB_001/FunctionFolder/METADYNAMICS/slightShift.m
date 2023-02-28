function Coords = slightShift(Coordinates,Lattice)

if length(Lattice) == 6
    Lattice = latConverter(Lattice);
end

% it is (a lot) easier to calculate the distances between ions, if we 
% first transform the Coordinates into a matrice, where each column represents
% the coordinates on one lattice vector (basis)
if ~isempty(find(size(Coordinates)==1))
    Coordinates_temp = Coordinates;
    Coordinates = zeros(3,size(Coordinates,1));
    Coordinates(:) = Coordinates_temp(:);
    Coordinates = Coordinates';
end

shiftDist = [-0.004,0,0]*inv(Lattice);

x_at = find(Coordinates(:,1)<0.1);
y_at = find(Coordinates(:,2)<0.1);
z_at = find(Coordinates(:,3)<0.1);
xy_at = find(Coordinates(:,1)<0.1  && Coordinates(:,2)<0.1);
xz_at = find(Coordinates(:,1)<0.1  && Coordinates(:,3)<0.1);
yz_at = find(Coordinates(:,2)<0.1  && Coordinates(:,3)<0.1);
xyz_at = find(Coordinates(:,1)<0.1 && Coordinates(:,2)<0.1 && Coordinates(:,3)<0.1);

C1  = Coordinates(x_at,:);
C1(:,1) = C1(:,1)+1;
C2 =Coordinates(y_at,:);
C2(:,2) = C2(:,2)+1;
C3 = Coordinates(z_at,:);
C3(:,3) = C3(:,3)+1;
C4 = Coordinates(xy_at,:);
C4(:,1) = C4(:,1)+1;
C4(:,2) = C4(:,2)+1;
C5 = Coordinates(xz_at,:);
C5(:,1) = C5(:,1)+1;
C5(:,3) = C5(:,3)+1;
C6 = Coordinates(yz_at,:);
C6(:,2) = C6(:,2)+1;
C6(:,3) = C6(:,3)+1;
C7 = Coordinates(xyz_at,:);
C7(:,1) = C7(:,1)+1;
C7(:,2) = C7(:,2)+1;
C7(:,3) = C7(:,3)+1;


COOOR = cat(1,Coordinates,C1,C2,C3,C4,C5,C6,C7);
flags = cat(1,ones(size(Coordinates,1),1),zeros(size(Coordinates(x_at,:),1)+size(Coordinates(y_at,:),1)+size(Coordinates(z_at,:),1)+size(Coordinates(xy_at,:),1)+size(Coordinates(xz_at,:),1)+size(Coordinates(yz_at,:),1)+size(Coordinates(xyz_at,:),1),1));

COOOR = COOOR*Lattice;

distOrigin = (sqrt(COOOR(:,1).^2 + COOOR(:,2).^2 + COOOR(:,3).^2));

[sDist,inds]=sort(distOrigin);

dosome= find(diff(sDist)<0.002);


doI = find(flags(inds(dosome)));

if ~isempty(doI)
Coordinates(inds(dosome(doI)),1) = Coordinates(inds(dosome(doI)),1)+shiftDist(1);
Coordinates(inds(dosome(doI)),2) = Coordinates(inds(dosome(doI)),2)+shiftDist(2);
Coordinates(inds(dosome(doI)),3) = Coordinates(inds(dosome(doI)),3)+shiftDist(3);
end

 Coords=Coordinates - floor(Coordinates);
  
