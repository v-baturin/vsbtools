function [lat, coord, numIons, chanAlist] = makeSurface(sur_lattice, sur_coordinate, sur_numIons, bulk_lattice, bulk_pos, bulk_numIons, vacuumSize)

  coord = [];
  chanAlist = [];
  if abs(det(bulk_lattice-sur_lattice))> 0.1
     disp('The lattice vector between bulk and surface is inconsistent');
     bulk_lattice
     sur_lattice
     quit
  elseif size(bulk_pos,1)~=sum(bulk_numIons)
     disp('numIons is inconsistent in substrate')
     bulk_numIons
     bulk_pos
     quit
  elseif size(sur_coordinate,1) ~= sum(sur_numIons)
     disp('numIons is inconsistent in surface')
     sur_numIons
     sur_coordinate
     quit
  else
     d1 = bulk_lattice(3,3);
     d2 = sur_lattice(3,3);
     d3 = vacuumSize;

     lat = bulk_lattice;
     lat(3,3) = d1+d2+d3;                                                     %total lattice

     if size(sur_coordinate,1) ~=0
        sur_coordinate(:,3) = (d1+sur_coordinate(:,3)*d2)/(d1+d2+d3); 
     end
   
     bulk_pos(:,3) = bulk_pos(:,3)*d1/(d1+d2+d3);
     numIons = sur_numIons+bulk_numIons;                                      %total Ions  
     count1 = 0;
     count2 = 0;
     for i = 1:length(numIons)
         coord = [coord; bulk_pos(count1+1:count1+bulk_numIons(i),:)];
         coord = [coord; sur_coordinate(count2+1:count2+sur_numIons(i),:)];   %total coords
         chanAlist = [chanAlist; zeros(bulk_numIons(i),1)];                   %total list
         chanAlist = [chanAlist; ones(sur_numIons(i),1)];
         count1 = count1+bulk_numIons(i);
         count2 = count2+sur_numIons(i);
     end
   end
