function extractStructures_PSO(max_num, weight)
% This function creates a list of structures with good fitness.
% max_num: The maximum number for output good structures
% weight:  for fingerprint analysis
% all the structure are loaded from PSO_STRUC

global PSO_STRUC

atomType = PSO_STRUC.atomType;
N = length(PSO_STRUC.POPULATION);
for i=1:N
    numIons    = PSO_STRUC.POPULATION(i).numIons;
    fitness(i) = PSO_STRUC.POPULATION(i).Fitness;
    enth(i)    = PSO_STRUC.POPULATION(i).Enthalpies(end)/sum(numIons);
end
[nothing, ranking] = sort(fitness);

% Remove all duplicates:
%--------------------------------------------------------------------------
% First step: fine filtering:
for i = 2 : N
    if PSO_STRUC.POPULATION(ranking(i)).ToCount > 0
        for j= 1 : i-1
            if PSO_STRUC.POPULATION(ranking(j)).ToCount > 0
                if SameStructure(ranking(i), ranking(j),PSO_STRUC)
                    PSO_STRUC.POPULATION(ranking(i)).ToCount = 0;
                    break;
                end
            end
        end
    end
end
%--------------------------------------------------------------------------
% Second step: rough filtering using A_order:
GoodList = zeros(max_num, 1); % Good Structure
item = 1;
GoodList(1) = 1; % Good Structure
for i = 2 : N
    if PSO_STRUC.POPULATION(ranking(i)).ToCount > 0
        same = 0;
        for j= 1 : i-1
            if PSO_STRUC.POPULATION(ranking(j)).ToCount > 0
                if SameStructure_order(ranking(i), ranking(j),PSO_STRUC)
                    PSO_STRUC.POPULATION(ranking(i)).ToCount = 0;
                    same = 1;
                    break;
                end
            end
        end
        if same == 0
            item = item + 1;
            GoodList(item) = i;
        end
    end
    if item == max_num
        break;
    end
end
%--------------------------------------------------------------------------

%safesave('../PSO.mat',PSO_STRUC); We don't save here

Good_num = item; %less than max_num
% output all accepted structures
fp1 = fopen('goodStructures', 'w');
fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes    fitness  SYMM\n');
fprintf(fp1,'                       (eV/atom)    (A^3/atom)    ()         \n');
unixCmd(['cat /dev/null > goodStructures_POSCARS']);  %Start to get new POSCAR
for i = 1 : Good_num
    j = ranking(GoodList(i));
    num    = PSO_STRUC.POPULATION(j).numIons;
    volume = PSO_STRUC.POPULATION(j).Vol/sum(num);
    symg   = PSO_STRUC.POPULATION(j).symg;
    lattice= PSO_STRUC.POPULATION(j).LATTICE;
    coor   = PSO_STRUC.POPULATION(j).COORDINATES;
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
        composition=[composition,blanks(shift(length(num)))];
    end
    
    fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %9.4f  %4d\n',...
        j, composition, enth(j), volume, fitness(j), symg);
    % POSCARS
    Write_POSCAR(atomType, j, symg, num, lattice, coor);
    unixCmd([' cat POSCAR      >> goodStructures_POSCARS']);
end
fclose(fp1);
