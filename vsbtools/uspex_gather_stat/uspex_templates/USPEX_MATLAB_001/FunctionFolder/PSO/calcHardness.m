function hardnesses = calcHardnessBondTypes(N_val, val1)

% calculates hardness for all current population structures
% R_val - covalent radii, N_val - number of valence electrons, val - valence
% we assume that R are given in A (thus rho = 0.37 in formula for nu)
% to include 'long' bonds we have to deal with graph connectivity (for 27=3x3x3 cells with cyclic boundary conditions)
% first we take into account all short bonds, then add long bonds between not connected subgraphs before graph is connected


global POP_STRUC
global ORG_STRUC

val = abs(val1);

L = length(POP_STRUC.POPULATION);

R_val = zeros(1,length(ORG_STRUC.atomType));

for i = 1 : length(ORG_STRUC.atomType)
    % R_val(i) = str2num(covalentRadius(ORG_STRUC.atomType(i)))
    s = covalentRadius(ceil(ORG_STRUC.atomType(i)));
    R_val(i) = str2num(s);
end


hardnesses = ones(1,L);
vect = zeros(1,3);
bonds_add = zeros(1,7);
same_bond = 0.05; % within this distance bond between same type of atoms are considered as bonds of the same time
max_bond = 5; % max possible distance for bonds, everything above that is ignored
small_bond = -0.37*log(ORG_STRUC.goodBonds); % max value that allows bond to be added when it doesn't increase connectivity, 0.05-0.1 seems reasonable
nL = 2; % number of 'layers' around the inner cell that we consider (should be either 1 or 2!, not more)


for p = 1 : L
    
    if POP_STRUC.POPULATION(p).FITNESSES(end) == 100000
        continue
    end
    
    lat = POP_STRUC.POPULATION(p).LATTICE;
    coords = POP_STRUC.POPULATION(p).COORDINATES;
    if length(lat) == 6
        lat = latConverter(lat);
    end
    N_i = sum(POP_STRUC.POPULATION(p).numIons);
    at_types = zeros(1,N_i);
    for k = 1 : N_i
        tmp = k;
        while tmp > 0
            at_types(k) = at_types(k) + 1;
            tmp = tmp - POP_STRUC.POPULATION(p).numIons(at_types(k));
        end
    end
    N_bonds1 = 0; % number of 'bonds'
    colors = zeros(2*nL+1, 2*nL+1, 2*nL+1, N_i); % used to check connectivity - all colors are the same in this case
    Ncolors = ((2*nL)^3)* N_i;
    for k1 = 1 : 2*nL+1
        for k2 = 1 : 2*nL+1
            for k3 = 1 : 2*nL+1
                for i = 1 : N_i
                    colors(k1,k2,k3,i) = i + ((2*nL+1)^2)*N_i*(k1-1) + (2*nL+1)*N_i*(k2-1) + N_i*(k3-1);
                end
            end
        end
    end
    colors1 = colors;
    it = 1;
    
    % check all different pairs of atoms
    % find all bonds not longer than 5A
    for i = 1 : N_i
        for j = i : N_i
            for k1 = -1 : 1
                for k2 = -1 : 1
                    for k3 = -1 : 1
                        if (i == j) & (k1 == 0) & (k2 == 0) & (k3 == 0)
                            continue;
                        end
                        vect(1) = coords(i,1) + k1 - coords(j,1);
                        vect(2) = coords(i,2) + k2 - coords(j,2);
                        vect(3) = coords(i,3) + k3 - coords(j,3);
                        delta = sqrt(sum((vect*lat).^2)) - R_val(at_types(i)) - R_val(at_types(j)); % wrong - sqrt(vect(1)^2*lat(1)^2 + vect(2)^2*lat(2)^2 + vect(3)^2*lat(3)^2) - R_val(a) - R_val(b);
                        if delta < max_bond
                            % add bond, keeping them sorted by distance
                            N_bonds1 = N_bonds1 + 1;
                            bonds_add(1) = i;
                            bonds_add(2) = j;
                            bonds_add(3) = delta;
                            bonds_add(4) = 0; % type
                            bonds_add(5) = k1;
                            bonds_add(6) = k2;
                            bonds_add(7) = k3;
                            
                            if N_bonds1 == 1
                                bonds = bonds_add;
                            else
                                bonds = vertcat(bonds, bonds_add);
                                for bb = 1 : N_bonds1-1
                                    if bonds(bb,3) > bonds_add(3)
                                        for bbb = N_bonds1 : -1 : bb+1
                                            bonds(bbb,:) = bonds(bbb-1,:);
                                        end
                                        bonds(bb,:) = bonds_add;
                                        break;
                                    end
                                end
                            end
                        end
                    end % k1
                end  % k2
            end   % k3
        end    % j
    end     % i
    
    % assign bond types
    bond_types = 0;
    N_bonds = 0;
    for bb = 1 : N_bonds1
        if bonds(bb,4) == 0;
            bond_types = bond_types + 1;
            bonds(bb,4) = bond_types;
            N_bonds = N_bonds + 1;
            if N_bonds == 1
                bonds_sorted = bonds(bb,:);
            else
                bonds_sorted = vertcat(bonds_sorted, bonds(bb,:));
            end
            for bbb = bb+1 : N_bonds1
                if bonds(bbb,3) > bonds(bb,3) + same_bond
                    break;
                end
                if (bonds(bbb,4) == 0) & (at_types(bonds(bbb,1)) == at_types(bonds(bb,1))) & (at_types(bonds(bbb,2)) == at_types(bonds(bb,2)))
                    bonds(bbb,4) = bond_types;
                    N_bonds = N_bonds + 1;
                    bonds_sorted = vertcat(bonds_sorted, bonds(bbb,:));
                end
            end
        end
    end
    
    % bonds_sorted
    % return;
    
    % add bonds by group, but only if this addition would increase the connectivity
    % important assumption - only neighbour cells are participating in connectivity
    % aka there is no shortest path connecting the atoms in the central unit cell that
    % goes more than 1 unit cell away form it
    bondTypesDensity = zeros(1,bond_types);
    h_tmp = ones(1,bond_types);
    it = 0;
    for bt = 1 : bond_types
        while (bonds_sorted(it+1,4) == bt) & (it < N_bonds)
            it = it + 1;
            a = at_types(bonds_sorted(it,1));
            b = at_types(bonds_sorted(it,2));
            bond_used = 0;
            for m1 = 1 : 2*nL+1
                for m2 = 1 : 2*nL+1
                    for m3 = 1 : 2*nL+1
                        n1 = m1 + bonds_sorted(it,5);
                        n2 = m2 + bonds_sorted(it,6);
                        n3 = m3 + bonds_sorted(it,7);
                        if (n1 > 2*nL+1) | (n1 < 1) | (n2 > 2*nL+1) | (n2 < 1) | (n3 > 2*nL+1) | (n3 < 1)
                            continue;
                        end
                        if (colors(nL+1,nL+1,nL+1,bonds_sorted(it,2)) == colors(nL+1+bonds_sorted(it,5), nL+1+bonds_sorted(it,6), nL+1+bonds_sorted(it,7),bonds_sorted(it,1))) & (bonds_sorted(it,3) > small_bond(a,b))
                            continue;
                        end
                        bond_used = 1;
                        cj = colors1(m1,m2,m3,bonds_sorted(it,2));
                        ci = colors1(n1,n2,n3,bonds_sorted(it,1));
                        for r1 = 1 : 2*nL+1
                            for r2 = 1 : 2*nL+1
                                for r3 = 1 : 2*nL+1
                                    for m = 1 : N_i
                                        if colors1(r1,r2,r3,m) == cj
                                            colors1(r1,r2,r3,m) = ci;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            bondTypesDensity(bt) = bondTypesDensity(bt) + bond_used;
            if bond_used == 0
                bonds_sorted(it,1) = 0;
                bonds_sorted(it,2) = 0;
            end
            if it == N_bonds
                break;
            end
        end
        % next step is needed since we don't apply color periodic boundary conditions
        for i = 1 : N_i
            for j = i+1 : N_i
                if colors1(nL+1,nL+1,nL+1,i) == colors1(nL+1,nL+1,nL+1,j)
                    for m1 = 1 : 2*nL+1
                        for m2 = 1 : 2*nL+1
                            for m3 = 1 : 2*nL+1
                                if colors1(m1,m2,m3,i) ~= colors1(m1,m2,m3,j)
                                    cj = colors1(m1,m2,m3,i);
                                    ci = colors1(m1,m2,m3,j);
                                    for r1 = 1 : 2*nL+1
                                        for r2 = 1 : 2*nL+1
                                            for r3 = 1 : 2*nL+1
                                                for m = 1 : N_i
                                                    if colors1(r1,r2,r3,m) == cj
                                                        colors1(r1,r2,r3,m) = ci;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        if it+1 < N_bonds
            if bonds_sorted(it+1,3) < max(max(small_bond))
                colors = colors1;
                continue;
            end
        end
        % check the connectivity, dvusvjaznostj dopustima
        colors = colors1;
        connected = 1;
        c1 = colors(1,1,1,1);
        c2 = 0;
        
        for k1 = nL : nL+2
            for k2 = nL : nL+2
                for k3 = nL : nL+2
                    for i = 1 : N_i
                        if (colors(k1,k2,k3,i) ~= c1) & (c2 == 0)
                            c2 = colors(k1,k2,k3,i);
                        elseif (colors(k1,k2,k3,i) ~= c1) & (colors(k1,k2,k3,i) ~= c2)
                            connected = 0;
                        end
                    end
                end
            end
        end
        if connected | it >= N_bonds
            for dumm = it+1 : N_bonds
                bonds_sorted(dumm,1) = 0;
                bonds_sorted(dumm,2) = 0;
            end
            break;
        end
    end
    
    % adjust coefficients nu so that their sum on the atom equals its N_val
    nu_factor = zeros(1,N_i);
    for k = 1 : N_i
        %  nu_full = sum(exp(-bonds((find((bonds(:,1) == k) | (bonds(:,2) == k)),3)/0.37)); % probably not compatible with oldMatlab on pegasus
        nu_full = 0;
        for m = 1 : N_bonds
            if (bonds_sorted(m,1) == k) | (bonds_sorted(m,2) == k)
                nu_full = nu_full + exp(-1*bonds_sorted(m,3)/0.37);
            end
        end
        nu_factor(k) = val(at_types(k))/nu_full;
        %  nu_factor(k) = 1;
    end
    CN = [6,3];
    it = 0;
    for bt = 1 : bond_types
        while (bonds_sorted(it+1,4) == bt) & (it < N_bonds)
            it = it + 1;
            if bonds_sorted(it,1) == 0
                if it == N_bonds
                    break;
                end
                continue;
            end
            a = at_types(bonds_sorted(it,1));
            b = at_types(bonds_sorted(it,2));
            delta = bonds_sorted(it,3);
            R_a = R_val(a) + delta/2;
            R_b = R_val(b) + delta/2;
            nu = exp(-delta/0.37);
            EN_a = 0.481*N_val(a)/R_a;
            EN_b = 0.481*N_val(b)/R_b;
            CN_a = val(a)/(nu*nu_factor(bonds_sorted(it,1)));
            CN_b = val(b)/(nu*nu_factor(bonds_sorted(it,2)));
            %   CN_a = bondTypesDensity(bt)/POP_STRUC.POPULATION(p).numIons(a);
            %   CN_b = bondTypesDensity(bt)/POP_STRUC.POPULATION(p).numIons(b);
            %   CN_a = CN(a);
            %   CN_b = CN(b);
            f_ab = 0.25*abs(EN_a - EN_b)/sqrt(EN_a*EN_b);
            X_ab = sqrt((EN_a*EN_b)/(CN_a*CN_b));
            h_tmp(bt) = h_tmp(bt)*X_ab*exp(-2.7*f_ab);
            if it == N_bonds
                break;
            end
        end
        if bondTypesDensity(bt) > 0
            h_tmp(bt) = bondTypesDensity(bt)*h_tmp(bt)^(1/bondTypesDensity(bt));
        end
    end
    for bt = 1 : bond_types
        hardnesses(p) = hardnesses(p)*h_tmp(bt);
    end
    differentTypes = sum(bondTypesDensity > 0);
    
    hardnesses(p) = 423.8*differentTypes*(hardnesses(p).^(1.0/differentTypes))/POP_STRUC.POPULATION(p).Vol(end) - 3.4;
end

%hardnesses = hardnesses;
