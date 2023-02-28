function MINHOP_extractStructures(max_num, weight)
% This function creates a list of structures with good fitness.
% max_num: The maximum number for output good structures
% weight:  for fingerprint analysis
% all the structure are loaded from USPEX_STRUC
% Last updated by Maxim Rakitin (2014/11/20)

global USPEX_STRUC

atomType = USPEX_STRUC.atomType;
N = length(USPEX_STRUC.POPULATION);
count = 1;
for i=1:N
    if USPEX_STRUC.POPULATION(i).ToCount > 0
        ID(count)   = i;
        numIons     = USPEX_STRUC.POPULATION(i).numIons;
        enth(count) = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(numIons);
        count = count + 1;
    end
end
[nothing, ranking] = sort(enth);

% Remove all duplicates:
%--------------------------------------------------------------------------
% First step: fine filtering:
if count >2
    for i = 2 : count-1
        if USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount > 0
            for j= 1 : i-1
                if USPEX_STRUC.POPULATION(ID(ranking(j))).ToCount > 0
                    if SameStructure(ID(ranking(i)), ID(ranking(j)),USPEX_STRUC)
                        USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount = 0;
                        break;
                    end
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
if count >2
    for i = 2 : count-1
        if USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount > 0
            same = 0;
            for j= 1 : i-1
                if USPEX_STRUC.POPULATION(ID(ranking(j))).ToCount > 0
                    if SameStructure_order(ID(ranking(i)), ID(ranking(j)),USPEX_STRUC)
                        USPEX_STRUC.POPULATION(ID(ranking(i))).ToCount = 0;
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
end
%--------------------------------------------------------------------------
%safesave('../USPEX.mat',USPEX_STRUC); We don't save here

Good_num = item; %less than max_num
% output all accepted structures
fp1 = fopen('goodStructures', 'w');
fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes     SYMM\n');
fprintf(fp1,'                       (eV/atom)    (A^3/atom)       \n');
unixCmd(['cat /dev/null > goodStructures_POSCARS']);  %Start to get new POSCAR
for i = 1 : Good_num
    j =  ID(ranking(GoodList(i)));
    num    = USPEX_STRUC.POPULATION(j).numIons;
    volume = det(USPEX_STRUC.POPULATION(j).LATTICE)/sum(num);
    symg   = USPEX_STRUC.POPULATION(j).symg;
    lattice= USPEX_STRUC.POPULATION(j).LATTICE;
    coor   = USPEX_STRUC.POPULATION(j).COORDINATES;
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
        composition=[composition,blanks(shift(length(num)))];
    end
    fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %4d\n',...
        j, composition, enth(ranking(GoodList(i))), volume, symg);
    % POSCARS
    Write_POSCAR(atomType, j, symg, num, lattice, coor);
    unixCmd([' cat POSCAR      >> goodStructures_POSCARS']);
    
end
fclose(fp1);
