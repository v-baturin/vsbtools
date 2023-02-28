function [ atom_position, atom_type, group_name, group_num,latt ] = Coordinate_read( filename )
%read the coordinate and atom label from cif files
% 

fid = fopen(filename,'r');
if fid == -1
    disp('file does not exist');
end 
str_label = [];
str_coor = [];
num=0;
num_atom = 0;
name1 = 0;
num1 = 0;
latt = zeros(6,1);
flag1 = 1;
flag2 = 0;
while ~feof(fid) 
    tmp = fgetl(fid);
    if findstr(tmp,'_cell_length_a')
       for i = 1 : 6
       latt(i) = sscanf(tmp, '%*s %g');
       tmp = fgetl(fid);
       end
    end
    
    
     while (~isempty(findstr(tmp,'_atom_site')) && isempty(findstr(tmp, 'atom_site_aniso'))) && flag1 ==1
            str_label = strvcat(str_label,tmp);
            num = num+1;
            tmp = fgetl(fid);
            if isempty((findstr(tmp,'_atom_site')))
                flag1 = 0;
                flag2 = 1;
                break;
            end
     end
    
    
     while flag2==1        
          if length(tmp) < 10 
              flag2 = 0;
             break;
          end              
           if  ~isempty(findstr(tmp,'loop')) || ~isempty(findstr(tmp,'#')) || isempty(tmp) || ~isempty(findstr(tmp,'atom_site'))
               flag2 = 0;
               break;
           end
              
           str_coor = strvcat(str_coor,tmp);
           num_atom = num_atom + 1;
           tmp = fgetl(fid);                         
     end
   
    if findstr(tmp, 'symmetry_space_group_name')
        group_name = sscanf(tmp, '%*s %s %s %s %s');
        name1 = name1 + 1;
    end
    
    if findstr(tmp, 'symmetry_Int_Tables_number')
        group_num = tmp(regexp(tmp,'\d'));
        group_num = str2num(group_num);
        num1 = num1 + 1;
    end
  
end

if name1 == 0          
    group_name = '0' % if there is no space group name, the default name is '0'
end
if num1 == 0
    group_num = 0;   % if there is no space group num, the default num is 0 
end

locate = 0;    % locate is the number of coordinate x in _atom_site 
element = 0;   % element is the number of element label in _atom_site 
num_label = size(str_label,1);
for i = 1 : size(str_label,1)
    if ~isempty (findstr(str_label(i,:),'_atom_site_fract_x'))
        locate = i;
    end
    if ~isempty (findstr(str_label(i,:),'_atom_site_type_symbol'))
        element = i;
    else if ~isempty (findstr(str_label(i,:),'_atom_site_label'))
            element = i;
        end
    end
end

% here we need to get the atom element symbol     

%tmp = 0;
atom_position = zeros(num_atom,3);   % for each atom, the position
atom_type = zeros(num_atom,1);
%element_name = zeros(num_atom,1);
for i = 1 : num_atom               % for each atom, read the atom type and the position
    position_1 = strtrim(str_coor(i,:));
    position_2 = regexp(position_1,'\s+', 'split');
    position_3 = char(position_2);
    element_name= getElementName(strtrim(position_3(element,:)));
    %tmp = strtrim(position_3(element,:));
    atom_type(i)= Atom_type(element_name);  % in position_3: 1 5 9 13 ... is the 
    atom_position(i,1)= getAtomCoords(position_3(locate,:));  % x: 2 6 10 14 ... 
    atom_position(i,2)= getAtomCoords(position_3(locate+1,:));  % y: 3 7 11 15 ...
    atom_position(i,3)= getAtomCoords(position_3(locate+2,:));   % z: 4 8 12 16 ...
end
atom_position = atom_position - floor(atom_position);

fclose(fid);
end




