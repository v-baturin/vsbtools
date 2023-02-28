function [a] = tensorProduct(t,b)
% tensor product in 3 dimensions
% a_ij = t_ijkl*b_kl

a = zeros(3,3);

for i = 1 : 3
 for j = 1 : 3
  for k = 1 : 3
   for l = 1 : 3
     a(i,j) = a(i,j) + t(i,j,k,l)*b(k,l);
   end
  end
 end
end
