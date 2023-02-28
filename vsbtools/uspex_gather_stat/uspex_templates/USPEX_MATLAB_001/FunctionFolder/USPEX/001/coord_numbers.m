function CN = coord_numbers(coord_v, numIons, lattice)

global ORG_STRUC

coord = coord_v*lattice;

%find type and covalent radius of each atom
for i = 1 : sum(numIons)
    for j = 1 : length(numIons)
        if i <= sum(numIons(1:j))
            typeAtom(i) = ORG_STRUC.atomType(j);
            R(i) = str2num(covalentRadius(typeAtom(i)));
            break;
        end
    end
end

%find coordiunation numbers of each atom
CN = [];
for i = 1 : sum(numIons)
   zn = [];
   for j = 1 : sum(numIons)
        if i==j
            zn(j)=0;
        else
            delta_r(j)=((coord(i,1)-coord(j,1))^2+(coord(i,2)-coord(j,2))^2+(coord(i,3)-coord(j,3))^2)^0.5;
            zn(j)=exp(-(delta_r(j) - R(i) - R(j))/0.23);
            if (typeAtom(i)==8) && (typeAtom(j)==8) %only for MOPAC Si-O clusters with MNDO
                zn(j)=exp(-((delta_r(j)+0.25) - R(i) - R(j))/0.23);
            end
        end
   end
   CN(i)=sum(zn)/max(zn);
end
