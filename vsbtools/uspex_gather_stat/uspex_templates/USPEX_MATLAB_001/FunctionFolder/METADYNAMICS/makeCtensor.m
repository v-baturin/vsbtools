function [S, C] = makeCtensor()

v = 0.26; % Poisson's ratio; threshold between brittle and elastic material

% http://ru.wikipedia.org/wiki/???????_??????
% makes Cijkl and Sijkl tensors from Cab:
% 1-v  v   v    0    0    0
%  v  1-v  v    0    0    0
%  v   v  1-v   0    0    0
%  0   0   0  0.5-v  0    0
%  0   0   0    0  0.5-v  0
%  0   0   0    0    0  0.5-v

% This Cab corresponds to isotropic material with Poisson's ratio = v 
% Young's modulus is choosen to 'normalize' the tensor for 'ideal'
% metadynamics case where pressure is in the one direction only and lattice is cubic

% rules for C: 11<=>1, 22<=>2, 33<=>3, 23,32<=>4, 13,31<=>5, 12,21<=>6
% for example, C_3122 = C0_52 
 rules = [1,6,5; 6,2,4; 5,4,3];
 C0 = [1-v,v,v,0,0,0; v,1-v,v,0,0,0; v,v,1-v,0,0,0; 0,0,0,0.5-v,0,0; 0,0,0,0,0.5-v,0; 0,0,0,0,0,0.5-v];
% C0 = C0/6;
 S0 = inv(C0);
 for i = 1 : 3
  for j = 4 : 6
   S0(i,j) = S0(i,j)/2;
   S0(j,i) = S0(j,i)/2;
  end
 end
 for i = 4 : 6
  for j = 4 : 6
   S0(i,j) = S0(i,j)/4;
  end
 end
% make Cijkl, Sijkl
 C = zeros(3,3,3,3);
 S = zeros(3,3,3,3);
 for i = 1 : 3
  for j = 1 : 3
   for k = 1 : 3
    for l = 1 : 3
     C(i,j,k,l) = C0(rules(i,j),rules(k,l));
     S(i,j,k,l) = S0(rules(i,j),rules(k,l));
    end
   end
  end
 end
 
 tmp = tensorProduct(S,[1 0 0; 0 0 0; 0 0 0]); %/norm([1 0 0 0 0 0]);
 S = S/tmp(1,1);

% makes Cijkl and Sijkl tensors from Cab:
% 2/3 1/3 1/3  0   0   0
% 1/3 2/3 1/3  0   0   0
% 1/3 1/3 2/3  0   0   0
%  0   0   0  1/6  0   0
%  0   0   0   0  1/6  0
%  0   0   0   0   0  1/6

%  C0 = [1,1,1,0,0,0; 0,1,1,0,0,0; 0,0,1,0,0,0; 0,0,0,1,0,0; 0,0,0,0,1,0;0,0,0,0,0,1];
% 0.25, since S0(4,4) = 4 s2323, etc
%  S0 = [1,-1,0,0,0,0; 0,1,-1,0,0,0; 0,0,1,0,0,0; 0,0,0,0.25,0,0; 0,0,0,0,0.25,0; 0,0,0,0,0,0.25];
% old version for C and S:
% 111000    1-1 0  0   0   0
% 011000    0 1-1  0   0   0
% 001000    0 0 1  0   0   0
% 000100    0 0 0 0.25 0   0
% 000010    0 0 0  0  0.25 0
% 000001    0 0 0  0   0  0.25
