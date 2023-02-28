function [coor, lat] = Read_MOPAC_Structure()

coor = [];
lat = [];

fp = fopen('input.arc','r');
while ~feof(fp)
   tmp = fgetl(fp);
   if ~isempty(findstr(tmp,'FINAL GEOMETRY OBTAINED'))
       tmp1 = fgetl(fp);
       tmp1 = fgetl(fp);
       tmp1 = fgetl(fp);
       tmp1 = fgetl(fp);
       while ~isempty(findstr(tmp1, '+1'))
          tmpin = sscanf(tmp1,'%*s %g %*s %g %*s %g %*s',[1,3]);
          if isempty(findstr(tmp1, 'Tv'))
             coor = [coor; tmpin];
          else
             lat = [lat; tmpin];
          end
          tmp1 = fgetl(fp);
       end
   end
end
fclose(fp);

if isempty(lat)
   coor = bsxfun(@minus, coor, mean(coor)); %Vectorized
   lat_len1 = max(coor(:,1)) - min(coor(:,1)) + 10;
   lat_len2 = max(coor(:,2)) - min(coor(:,2)) + 10;
   lat_len3 = max(coor(:,3)) - min(coor(:,3)) + 10;
   lat = diag([lat_len1, lat_len2, lat_len3]);
   coor = bsxfun(@plus, coor, [lat_len1, lat_len2, lat_len3]/2);
elseif size(lat,1) ~=3
   disp('LATTICE is wrong......')
   quit
end
coor = coor/lat;
coor = coor - floor(coor);
