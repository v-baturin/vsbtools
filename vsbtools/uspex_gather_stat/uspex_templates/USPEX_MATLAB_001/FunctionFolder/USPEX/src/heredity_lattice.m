function Lattice = heredity_lattice(lat1, lat2, frac, constLat)
if ~constLat
    fracLattice = rand(1);
    lat1 = latConverter(latConverter(lat1));
    lat2 = latConverter(latConverter(lat2));
    temp_potLat = fracLattice*lat1 + (1-fracLattice)*lat2;
    volLat = det(temp_potLat);
    if volLat < 0
        temp_potLat = -1*temp_potLat;
    end

    latVol = det(lat1)*frac + det(lat2)*(1-frac);
    ratio = latVol/volLat;
    temp_potLat = latConverter(temp_potLat);
    temp_potLat(1:3)= temp_potLat(1:3)*(ratio)^(1/3);
    Lattice = latConverter(temp_potLat);
else
    Lattice = (lat1+lat2)/2;
end

