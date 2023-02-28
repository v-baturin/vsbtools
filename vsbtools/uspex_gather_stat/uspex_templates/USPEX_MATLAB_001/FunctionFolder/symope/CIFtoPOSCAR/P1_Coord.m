function [ position, N_1] = P1_Coord( atom_position, atom_type,N_type)
% For the P1 Space group, we output the coordinates of the atoms
% position is the coordinates of the atom according to the atom order in cif 
% N_1 is that the number of atoms for each type 

Natom = size(atom_position,1);
position_order = [];
N_1 = zeros(size(N_type,1),1);   % how may atoms for each type  Ex: 1 2 12 6 ( 1 Ca atoms, 2 Cl atoms, 12 H atoms, 6 O atoms)

for i = 1 : size(N_type,1)
    for j = 1 : Natom
        if atom_type(j) == N_type(i)
           N_1(i,1) = N_1(i,1) + 1 ;
           position_order(i,N_1(i,1),:)=atom_position(j,:);
        end
    end
end

position = zeros(Natom,3);
tmp = 1;
for j = 1 : size(N_type,1)
    for k = 1 : N_1(j,1)
        position(tmp,:) = position_order(j,k,:);  % according to the order of the N_type, output every atom's coordinate 
        tmp = tmp + 1;
    end
end


end

