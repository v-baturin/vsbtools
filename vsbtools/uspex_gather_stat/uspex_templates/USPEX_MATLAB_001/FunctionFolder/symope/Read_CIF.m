function [lat, coor, operation, Error] = Read_CIF(filename)
lat       = [];
coor      = [];
operation = [];
Error     = 0;

fp = fopen(filename,'r');
while ~feof(fp)
  tmp = fgetl(fp); % repeats input file
  if findstr(tmp, 'bombed')
     Error = 1;
     break;

% _cell_length_a      15.20090
% _cell_length_b      10.31400
% _cell_length_c       8.85900
% _cell_angle_alpha   90.00000
% _cell_angle_beta    85.86133
% _cell_angle_gamma   90.00000
  elseif  ~isempty(findstr(tmp,'_cell_'))
      tmp = regexp( tmp, '\s+', 'split');
      lat = [lat; str2num(tmp{2})];

%_space_group_symop_operation_xyz
%x,y,z
%-x,y+1/2,-z+1/2
%-x,-y,-z
%x,-y+1/2,z+1/2
  elseif ~isempty(findstr(tmp,'symop_operation_xyz'))
      tmp = fgetl(fp); % repeats input file
      while ~isempty(findstr(tmp, 'x'))
           operation = [operation; GetMatrix(tmp)];
           tmp = fgetl(fp); % repeats input file
      end

%_atom_site_occupancy
%A1 A  -0.19204   0.46823  -0.08784   1.00000
  elseif ~isempty(findstr(tmp,'_atom_site_occupancy'))
      tmp = fgetl(fp); % repeats input file
      while ~isempty(findstr(tmp, '1.0000'))
           tmp = regexp( tmp, '\s+', 'split');
           xyz = [str2num(tmp{3}), str2num(tmp{4}), str2num(tmp{5})];
           coor = [coor; xyz];
           tmp = fgetl(fp); % repeats input file
      end
  end
end

fclose(fp);

if isempty(lat) || isempty(coor)
   disp([pwd ': ' filename ' is broken']);
   Error = 1;
   %quit
else
   lat(4:6)=lat(4:6)*pi/180;
end

