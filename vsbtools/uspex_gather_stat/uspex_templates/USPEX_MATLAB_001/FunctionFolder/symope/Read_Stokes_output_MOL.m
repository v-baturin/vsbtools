function [Operation, Num, candidate, lattice, errorMsg] = Read_Stokes_output_MOL(filename)

%--This code can be used to read 2D/3D crystals

lattice = [1 1 1 0 0 0]';
candidate = [];
errorMsg  = 0;
nsym      = [];
numIons   = [];
Num       = [];
Operation = [];
if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end

fp = fopen(filename,'r');
if OctaveMode %Octave 3.4 does not support regexp('split')
   while ~feof(fp)
       tmp = fgetl(fp); % repeats input file
       if ~isempty(findstr(tmp,'Space-group symmetry'))
           tmp = strsplit( tmp, ':');
           nsym = tmp{2};
       elseif ~isempty(findstr(tmp,'atoms of each type: '))
           tmp = strsplit( tmp, ':');
           numIons = str2num(tmp{2});
       elseif ~isempty(findstr(tmp,'operation '))  %must space
           tmp = strrep( tmp, 'operation', '');
           Operation = [Operation; str2num(tmp)];
       elseif ~isempty(findstr(tmp,'NUMBER'))
           Num = [Num; str2num(fgetl(fp))];
       elseif ~isempty(findstr(tmp,'cell parameters'))
           tmp = strsplit( tmp, ':');
           lat = str2num(tmp{2});
           lattice = lat';
       elseif ~isempty( findstr(tmp, 'atomic parameters') )
           tmp = fgetl(fp);
           while ~feof(fp)
                tmp = fgetl(fp);
                if ~isempty(tmp)
                   tmp =strsplit( tmp, '(');
                   tmp = str2num(tmp{1});
                   candidate = [candidate; tmp(3:5)];
                end
           end
       elseif ~isempty( findstr(tmp, 'error') )
           if ~isempty( findstr(tmp, 'not possible') )
              errorMsg = 2;  %mismatch between N_atom and space group
           else
              errorMsg = 1;  %general errorMsg
           end
           %disp(WrongMes);
           break;
       end
   end

else %MATLAB

   while ~feof(fp)
       tmp = fgetl(fp); % repeats input file
       if ~isempty(findstr(tmp,'Space-group symmetry'))
           tmp = regexp( tmp, ':', 'split');
           nsym = tmp{2};
       elseif ~isempty( findstr(tmp, 'error') )
           if ~isempty( findstr(tmp, 'not possible') )
              errorMsg = 2;  %mismatch between N_atom and space group
           else
              errorMsg = 1;  %general errorMsg
           end
           break;
       elseif ~isempty(findstr(tmp,'atoms of each type: '))
           tmp = regexp( tmp, ':', 'split');
           numIons = str2num(tmp{2});
       elseif ~isempty(findstr(tmp,'operation '))
           tmp = regexp( tmp, '\s+', 'split');
           tmp = [tmp{3} ' ' tmp{4} ' ' tmp{5}];
           Operation = [Operation; str2num(tmp)];
       elseif ~isempty(findstr(tmp,'NUMBER'))
           Num = [Num; str2num(fgetl(fp))];
       elseif ~isempty(findstr(tmp,'cell parameters'))
           tmp = regexp( tmp, ':', 'split');
           lat = str2num(tmp{2});
           lattice = lat';
       elseif ~isempty( findstr(tmp, 'atomic parameters') )
           tmp = fgetl(fp);
           while ~feof(fp)
                tmp = fgetl(fp);
                if ~isempty(tmp)
                   tmp = regexp( tmp, '\s+', 'split');
                   tmp = [tmp{4} ' ' tmp{5} ' ' tmp{6}];
                   candidate = [candidate; str2num(tmp)];
                end
           end
       end
   end

end

fclose(fp);

if errorMsg == 0
   if sum(numIons) ~= size(candidate,1)
      disp(['Stokes output error: atomic coordinates are insufficient']);
   else
      if size(Num,2)>1
         Num = (sum(Num'))';
      end
      Num_tmp =[0;Num(1:length(Num)-1)];
      Num = Num - Num_tmp;
   end
end
