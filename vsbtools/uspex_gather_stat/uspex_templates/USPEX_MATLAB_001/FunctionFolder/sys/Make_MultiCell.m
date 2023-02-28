function [S_coor, S_lat, S_numIons, trans] = Make_MultiCell(coor, lat, numIons, vector)

%To make the supercell according to the vector matrix 
%[a1,a2;a3,a4]
coor = coor - floor(coor);  %exactly[0,1];
Grid = zeros(4,2);
Grid(2,:) = [vector(1),vector(2)];
Grid(3,:) = [vector(1)+vector(3),vector(2)+vector(4)];
Grid(4,:) = [vector(3),vector(4)];
trans = [vector(1) vector(2) 0;vector(3) vector(4) 0; 0 0 1];

X_Max = max(Grid(:,1));
Y_Max = max(Grid(:,2));
X_Min = min(Grid(:,1));
Y_Min = min(Grid(:,2));
Matrix = SuperMatrix(X_Min, X_Max-1, Y_Min, Y_Max-1, 0, 0);
N_Size = (X_Max-X_Min)*(Y_Max-Y_Min);
N_atom = size(coor,1);
tmp      = repmat(coor, [1, N_Size]);
%  ATOM: a1 b1 c1   ---->    a1 b1 c1
%        a2 b2 c2   ---->    a2 b2 c2
%                            a1 b1 c1
%                            a2 b2 c2
tmp_coor = reshape(tmp',[3,N_atom*N_Size])';
S_Matrix = repmat(Matrix, [N_atom,1]);
tmp_coor = (tmp_coor + S_Matrix)*lat;
S_lat    = trans*lat;
tmp_coor = tmp_coor/S_lat;
tmp_coor = tmp_coor - floor(tmp_coor);

tmp_numIons = size(tmp_coor,1);
%disp(['MakeMulti: numIons: ' num2str(tmp_numIons)])
[S_coor, tmp_numIons] = RemoveDup(tmp_coor, tmp_numIons, 0.01);
%disp(['MakeMulti: numIons: ' num2str(tmp_numIons)])

S_numIons = numIons*det(trans);
S_lat = latConverter(latConverter(S_lat));
if det(trans) ~= size(S_coor,1)/N_atom
   N_atom
   disp([num2str(det(trans)), '====' num2str(size(S_coor,1)/N_atom)]);
   disp('1.0000')
   disp(num2str(S_lat))
   disp('Al O')
   disp(num2str(S_numIons))
   disp('Direct')
   disp(num2str(S_coor))
   vector
   lat
   coor
   quit
end
