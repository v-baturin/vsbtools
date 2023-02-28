function Ope = GetMatrix(string)
Ope = zeros(4,3);
tmp = regexp(string, ',', 'split');
ref = {'x','y','z'};
for i = 1:3
   for j = 1:3
       pos = findstr(tmp{i},ref{j});
       if ~isempty(pos)
          Ope(i,j) = 1;
          tmp{i}(pos)=[];
          if (pos>1)
             if (tmp{i}(pos-1)=='-')
                Ope(i,j) = -1;
             end
             tmp{i}(pos-1)=[];
          end
       end
   end
   if ~isempty(tmp{i})
      Ope(4,i) = str2num(tmp{i});
   end
end

