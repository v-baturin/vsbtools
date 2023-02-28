function lat = perturbCell(lattice)

global ORG_STRUC

% perturbs lattice to break the symmetry
lat = lattice; 

if size(ORG_STRUC.lattice,1) ~= 3   % Cell hasn't been specified by user
    lat_6 = latConverter(lattice);
    ABC = 0.02*lat_6(1:3).*(rand(3,1)-0.5);
    angle = (rand(3,1)-0.5)*pi/90;
    lat_6 = lat_6 + [ABC; angle];
    lat = latConverter(lat_6);
end

