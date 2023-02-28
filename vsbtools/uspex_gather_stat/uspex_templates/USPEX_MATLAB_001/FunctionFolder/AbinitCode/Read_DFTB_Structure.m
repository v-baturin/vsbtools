function [COORDS, LATT] = Read_DFTB_Structure()
%This rountine is to read crystal structure 
%File: geo_end.gen 

fp = fopen('geo_end.gen');
tmp = fgetl(fp); 
natom = str2num(tmp(1:6));
tmp = fgetl(fp); 
COORDS = zeros(natom, 3);
LATT   = zeros(3, 3);
for i = 1:natom
   tmp = str2num(fgetl(fp)); 
   COORDS(i,:) = tmp(3:5);
end
tmp = fgetl(fp); 
for i = 1:3
   tmp = fgetl(fp); % system description
   LATT(i,:) = str2num(tmp);
end
fclose(fp);

