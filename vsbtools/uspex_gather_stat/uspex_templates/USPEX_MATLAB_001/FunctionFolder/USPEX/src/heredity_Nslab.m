function Nslabs = heredity_Nslab(lat, numIons, corr, dimension)

% determine the thickness of the cell when cutting in a given direction 
if dimension == 1
    L = abs(det(lat)/norm(cross(lat(2,:),lat(3,:))));
elseif dimension == 2
    L = abs(det(lat)/norm(cross(lat(1,:),lat(3,:))));
else
    L = abs(det(lat)/norm(cross(lat(1,:),lat(2,:))));
end
% characteristic length ~= = 0.5*(V/N)^1/3
Lchar = 0.5*power(abs(det(lat))/sum(numIons), 1/3);
Nslabs = round(L/(Lchar+(L-Lchar)*(cos(corr*pi/2))^2));
%disp(['Number of slabs = ' num2str(Nslabs1) ' ' num2str(Nslabs2)])

