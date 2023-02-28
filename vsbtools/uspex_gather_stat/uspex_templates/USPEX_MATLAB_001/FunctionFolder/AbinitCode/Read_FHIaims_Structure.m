function [coor, lat] = Read_FHIaims_Structure()
%% In FHI-081213 geometry.in.next_step automatically will be created but
%% for FHI-081219 user need to specify restart_relaxations .true.

%% This condition was applied because sometimes which systems is small
%% USPEX creates very good structures which are the same with the relaxed one
%% and FHI finishes without changing the relaxed structure, thus  
%% geometry.in.next_step won't be created.
if ~exist('geometry.in.next_step') 
    unixCmd('cp geometry.in  geometry.in.next_step');  
end
%----------------------------------------------------------------------

command = ['grep "lattice_vector" geometry.in.next_step | cut -c16-62 '];
[nothing, lattice] = unix(command);

lat = [];

if ~isempty(lattice)
   lat = str2num(lattice);
end

command = ['grep "atom" geometry.in.next_step | cut -c6-53 '];
[nothing, coordinates] = unix(command);
coor = str2num(coordinates);

if isempty(lat)
   coor = bsxfun(@minus, coor, mean(coor)); %Vectorized
   lat_len1 = max(coor(:,1)) - min(coor(:,1)) + 10;
   lat_len2 = max(coor(:,2)) - min(coor(:,2)) + 10;
   lat_len3 = max(coor(:,3)) - min(coor(:,3)) + 10;
   lat = diag([lat_len1, lat_len2, lat_len3]);
   coor = bsxfun(@plus, coor, [lat_len1, lat_len2, lat_len3]/2);
end

coor = coor/lat;
coor = coor - floor(coor);
