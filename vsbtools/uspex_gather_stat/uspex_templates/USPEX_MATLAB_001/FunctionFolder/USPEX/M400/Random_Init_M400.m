function potentialOffspring = Random_Init_M400(angles_num)
% $Rev$
% $Author$
% $Date$

sec_parts = combineSecStructs(angles_num);
sequence  = randperm(7);

disp(' ');

potentialOffspring = [];
for i=1:size(sequence, 2)
    for j=1:sec_parts(sequence(i))
        [phi, psi, name] = secStructs(sequence(i));
        disp(name);
        potentialOffspring = [potentialOffspring; [double(phi), double(psi)]];
    end
end

end
