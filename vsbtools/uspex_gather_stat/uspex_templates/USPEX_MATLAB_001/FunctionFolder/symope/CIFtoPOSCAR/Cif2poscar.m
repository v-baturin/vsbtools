function Natoms = Cif2poscar()
%This function is to change the cif to POSCAR.
% 

cif_file    = 'symmetrized.cif';
poscar_file = 'outputPOSCAR';

fid = fopen(cif_file,'r');   % Open cif file
if fid == -1
    disp('file does not exist');
    %return
end

fp = fopen(poscar_file, 'w');


% here we begin to read the types and coordiantes of atoms, the space group, space group num and lattice in cif 
[atom_position, atom_type, group_name, group_num,lattice] = Coordinate_read(cif_file);
if ~isempty(atom_position)
% These angles in degrees are for a POSCAR header:
lattice_deg(1) = lattice(4);           % alpha
lattice_deg(2) = lattice(5);           % beta
lattice_deg(3) = lattice(6);           % gamma

% [MR]: convert angles from degrees to radians since MATLAB works with radians:
lattice(4) = lattice(4)*pi/180;        % alpha
lattice(5) = lattice(5)*pi/180;        % beta
lattice(6) = lattice(6)*pi/180;        % gamma

Lattice = latConverter(lattice');

% Lattice parameter is written to POSCAR
fprintf(fp, 'EA0 %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f \n', lattice(1:3), lattice_deg(1:3)); 
fprintf(fp, '%6.4f \n', 1.0);
for latticeLoop = 1 : size(Lattice,1)
  fprintf(fp, '%12.6f %12.6f %12.6f\n', Lattice(latticeLoop,:));
end



% delete the repeated atoms in cif
atom_position_np = atom_position(1,:);
atom_type_np = atom_type(1,:);
ntmp = 0;
npPOS = 1;
for i = 2 : size(atom_position,1)   
    for j = 1 : npPOS
        if  norm((atom_position(i,:)-atom_position_np(j,:))*Lattice) > 0.001
            ntmp = ntmp + 1;
        end
    end
    if ntmp == npPOS
       npPOS = npPOS+1;
       atom_position_np = [atom_position_np;atom_position(i,:)];
       atom_type_np = [atom_type_np;atom_type(i,:)];
    end 
    ntmp = 0;
end

% how many types of atom in the file using the order in cif. Ex: 20 17 1 8 (Ca Cl H O)
N_type_1 = length (unique(atom_type_np));
N_type = zeros(N_type_1,1); 
N_type(1) = atom_type_np(1);
num = 1;
num_tmp = 0;
for i = 2 : size(atom_type_np,1)
    for j = 1 : num
       if atom_type_np(i) ~= N_type(j)
          num_tmp = num_tmp +1 ;
      end
    end
    if num_tmp == num 
        num = num +1;
        N_type(num) = atom_type_np(i,1);
    end
    num_tmp = 0;   
end

N_atom =size(atom_position_np,1); % how many norepeated atoms in cif file 
 


% For the space group P1, we get the position and N_1
flag1 = 1; % evaluate whether the symmetry is P1
if strcmp(group_name,'''P1''')
   % (group_name =='P1')
    %|| (group_num == 1)
    flag1 = 0;
    [position,N_1] = P1_Coord(atom_position_np,atom_type_np, N_type); % for P1, we get the atom position and N_1.   N_1: the number of atoms in each type of atom
    
end
 
% evaluate whether there are the transformation operator
flag2 = 1;
flag3 =1;
while(flag1)
   transStr = callAWK('transform.awk', cif_file, '');
   % [nothing, transStr] = unix(['awk -f transform.awk ' cif_file]);
    if isempty(transStr)
        flag2 = 0; % flag2 : evaluate whether there is a transformation operator 
        disp('There is no transformation information in cif file') 
        flag3 = 0; % if there is no operator, flag3 will control the writing of the POSCAR and it will not continue to write the coordinates in the POSCAR 
    end


% For space group ~= 1 and there are space group operators 
while (flag2)
% here we begin to read the transformation operator in cif
trans_1 = strtrim(transStr);% delet the empty space at the beginning and the end
if (~isletter(trans_1(1)) && isspace(trans_1(2))) || trans_1(1) == ''''
    trans_1(1)= ' ';
end
for i = 2 :length(trans_1)              % delete quotes and number in the transformation operator
   if ~isletter(trans_1(i)) && isspace(trans_1(i-1)) && isspace(trans_1(i+1))
       trans_1(i)= ' ';
   else if trans_1(i) == ''''
           trans_1(i) = ' ';
       end
   end
end

trans_1 = strtrim(trans_1);

% here we begin to transform the operator to matrix  ex:  x y z is the operator. [1 0 0;0 1 0; 0 0 1] is the matrix 
trans_2 = regexp(trans_1, '\s+', 'split');  % split the array accoring to the empty space and trans_2 is a cell
N_trans =length(trans_2)/3;
trans_Matr = vpa(zeros(length(trans_2),3));% trans_Matr: the transformation matrix: a11, a12, a13,
trans_tran = vpa(zeros(length(trans_2),1));  % trans_tran: the translation matrix:
%trans_tran = zeros(length(trans_2),1);  % trans_tran: the translation matrix:
trans_num  = vpa(zeros(length(trans_2),1));  % the number of the transformation matrix
syms x y z m n;
trans_3 = trans_2;
%trans_3 = sym(trans_3);
%trans_3 = sym(trans_2);
%m = trans_3;
for i = 1 :length(trans_2) 
    tmp = coef(trans_3(i),'xyz');
%    tmp = coef(sym(m(i)),'xyz');
    trans_Matr(i,1)=tmp(1);  % trans_Matr is all the transformation matrix
    trans_Matr(i,2)=tmp(2);
    trans_Matr(i,3)=tmp(3);
    const = trans_3(i)-tmp(1)*x-tmp(2)*y-tmp(3)*z;
%    const = m(i)-sym(tmp(1))*x-sym(tmp(2))*y-sym(tmp(3))*z;
%    trans_tran(i) = const;  % trans_tran is the translation matrix ex [0.5;0.5;0] 
    trans_tran(i) = vpa(const);  % trans_tran is the translation matrix ex [0.5;0.5;0] 
    trans_num(i) = ceil(i/3);
end 
 

  
% here we begin to use the transformation matrix to check all the atoms in order to evaluate whether there is a new atom generated

New_position_1 = atom_position_np;
New_atomtype_1 = atom_type_np;
newPOS2 = N_atom;
ntmp2=0;
for i = 1 : N_atom
    for j = 1 : N_trans
        Matrix_1 = trans_Matr(3*j-2:3*j,:);
        Ini = atom_position_np(i,:);
        tmp = Matrix_1*Ini' + trans_tran(3*j-2:3*j,:);
        tmp = tmp';
        tmp = tmp - floor(tmp);
        %if ( (abs(tmp(1)-Ini(1)) < 0.0001) || (abs(tmp(1)+Ini(1)-1)< 0.0001) ) &&  ( (abs(tmp(2)-Ini(2)) < 0.0001) || (abs(tmp(2)+Ini(2)-1)< 0.0001) ) && ( (abs(tmp(3)-Ini(3)) < 0.0001) || (abs(tmp(3)+Ini(3)-1)< 0.0001) )
        for k = 1 : newPOS2
            if ( norm((tmp-New_position_1(k,:))*Lattice) > 0.001 )
                    ntmp2 = ntmp2 + 1;
           end
        end
        if ntmp2 == newPOS2
           newPOS2 = newPOS2 + 1;
           New_position_1 = [New_position_1 ; tmp];
           New_atomtype_1 = [New_atomtype_1; atom_type_np(i)];
        end

        ntmp2 = 0;
    end
end


N_1 = zeros(size(N_type,1),1);   % how may atoms for each type  Ex: 1 2 12 6 ( 1 Ca atoms, 2 Cl atoms, 12 H atoms, 6 O atoms)

% here output the atom coordiante and type with the same sequence of N_type
all_position = New_position_1;
all_atomtype = New_atomtype_1;
all_Natom =size(all_atomtype,1);    
position_order = [];
for i = 1 : size(N_type,1)
    for j = 1 : all_Natom
        if all_atomtype(j) == N_type(i)
           N_1(i,1) = N_1(i,1) + 1 ;
           position_order(i,N_1(i,1),:)=all_position(j,:);
        end
    end
end

position = zeros(all_Natom,3);
tmp = 1;
for j = 1 : size(N_type,1)
    for k = 1 : N_1(j,1)
        position(tmp,:) = position_order(j,k,:);  % according to the order of the N_type, output every atom's coordinate 
        tmp = tmp + 1;
    end
end


flag2 = 0;

end

flag1 = 0;

end

while(flag3)
for atomtypeLoop = 1 : size(N_type,1)
   if atomtypeLoop == size(N_type,1)
      fprintf(fp, '%4s\n', megaDoof(N_type(atomtypeLoop)));
   else
      fprintf(fp, '%4s', megaDoof(N_type(atomtypeLoop))); 
   end
end



for atomnumLoop = 1 : size(N_1,1)
   if atomnumLoop == size(N_1,1)
      fprintf(fp, '%4d\n', N_1(atomnumLoop,:));
   else
      fprintf(fp, '%4d', N_1(atomnumLoop,:)); 
   end
end


fprintf(fp, '%6s\n', 'Direct');

for coordLoop = 1 : size(position,1)
  fprintf(fp, '%12.6f %12.6f %12.6f\n', position(coordLoop,:));
end



flag3 = 0;
end
Natoms=length(atom_position);
else 
Natoms=0;
end


fclose(fid);
fclose(fp);

