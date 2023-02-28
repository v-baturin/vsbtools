function [ atomic_number ] = Atom_type( atomsymbol )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

  for j = 1 : 105
        if  strcmp(lower(atomsymbol), lower(megaDoof(j)))
             atomic_number = j;
             break;
        end
  end
    
end

