function [coor, lat] = Read_Gaussian_Structure()
global ORG_STRUC
[nothing, line] = unix('grep -n "Standard orientation:" Gauss_output | tail -1 | awk ''{print $1}''');
line = str2num(line(1:(length(line)-2)));
NUM = sum(ORG_STRUC.numIons(:));
unix(['sed -n ' num2str(line+5) ',' num2str(line+5+NUM-1) 'p Gauss_output > OUT']);
[nothing, coordinate] = unix('awk ''{print $4, $5, $6}'' OUT');
coor = str2num(coordinate);
lat_len = 30;
lat = [lat_len 0 0; 0 lat_len 0; 0 0 lat_len];
coor = bsxfun(@minus, coor, mean(coor)); %Vectorized
coor = coor + lat_len/2;
coor = coor/lat;
coor = coor - floor(coor);
unix('rm OUT');
