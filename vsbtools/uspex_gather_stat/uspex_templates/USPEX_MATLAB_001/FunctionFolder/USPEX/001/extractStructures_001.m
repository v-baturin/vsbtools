function extractStructures(max_num)
% This function creates a list of structures with good fitness.
% max_num: The maximum number for output good structures
% all the structure are loaded from USPEX_STRUC
% Last updated by Qiang Zhu (2014/02/18)

global USPEX_STRUC
global ORG_STRUC
global POP_STRUC

%In principle, weigth can be also obtained from USPEX, will do it later
atomType = USPEX_STRUC.SYSTEM(1).atomType;
N = length(USPEX_STRUC.POPULATION);

if ~exist('__goodStructures')
    mkdir('__goodStructures')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------%
%------------------- case of two-component (binary) clusters ---------------------%
%---------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(ORG_STRUC.numIons,2) == 2

for comp_x = ORG_STRUC.numIons(1,1) : ORG_STRUC.numIons(2,1)
for comp_y = ORG_STRUC.numIons(1,2) : ORG_STRUC.numIons(2,2)
  
  numions = [comp_x,comp_y];
  fitness = [];
  enth = [];
  number = []; % numner of current structure in USPEX_STRUC (array)
  nn = 0;      % the total number of structures in USPEX_STRUC with current composition
  for i=1:N
    if (USPEX_STRUC.POPULATION(i).numIons == numions)
      nn = nn+1;
      number(nn) = i;
      fitness(nn) = USPEX_STRUC.POPULATION(i).Enthalpies(end);
      enth(nn)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(numions);
    end
  end
  [nothing, ranking] = sort(fitness);

  % Remove all duplicates:
  %--------------------------------------------------------------------------
  GoodList = zeros(max_num, 1); % Good Structure
  item = 1;
  GoodList(1) = 1; % Good Structure
  for i = 2 : nn
    ii = number(ranking(i)); %number in USPEX_STRUC
    if USPEX_STRUC.POPULATION(ii).ToCount > 0
        same = 0;
        for j= 1 : i-1
            jj = number(ranking(j)); %number in USPEX_STRUC
            if USPEX_STRUC.POPULATION(jj).ToCount > 0
                if SameStructure(ii,jj,USPEX_STRUC)
                    USPEX_STRUC.POPULATION(ii).ToCount = 0;
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
  %safesave('../USPEX.mat',USPEX_STRUC); We don't save here

  Good_num = item; %less than max_num
  % output all accepted structures
  string_comp = ['composition_' num2str(comp_x) '_' num2str(comp_y)];
  fpath  = ['__goodStructures/' string_comp];
  fp1 = fopen(fpath, 'w');
  unixCmd(['cat /dev/null > ' fpath]);  %Start to get new POSCAR
  unixCmd(['cat /dev/null > ' fpath '_POSCARS']);
  %fp2 = fopen('goodStructures_POSCARS', 'w');
  fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes   Enthalpies SYMM\n');
  fprintf(fp1,'                       (eV/atom)    (A^3/atom)    (eV)        \n');
  for i = 1 : Good_num
    jj = number(ranking(GoodList(i)));
    num    = USPEX_STRUC.POPULATION(jj).numIons;
    volume = USPEX_STRUC.POPULATION(jj).Vol/sum(num);
    symg   = USPEX_STRUC.POPULATION(jj).symg;
    lattice= USPEX_STRUC.POPULATION(jj).LATTICE;
    coor   = USPEX_STRUC.POPULATION(jj).COORDINATES;
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
        composition=[composition,blanks(shift(length(num)))];
    end
    
    j = ranking(GoodList(i));
    fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %9.4f  %4d\n',...
        jj, composition, enth(j), volume, fitness(j), symg);
    % POSCARS
    Write_POSCAR(atomType, jj, symg, num, lattice, coor);
    unixCmd([' cat POSCAR      >> ' fpath '_POSCARS']);  
  end
  fclose(fp1);

end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------%
%------------------- case of three-component (binary) clusters ---------------------%
%---------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(ORG_STRUC.numIons,2) == 3

for comp_x = ORG_STRUC.numIons(1,1) : ORG_STRUC.numIons(2,1)
for comp_y = ORG_STRUC.numIons(1,2) : ORG_STRUC.numIons(2,2)
for comp_z = ORG_STRUC.numIons(1,3) : ORG_STRUC.numIons(2,3)
  
  numions = [comp_x,comp_y, comp_z];
  fitness = [];
  enth = [];
  number = []; % numner of current structure in USPEX_STRUC (array)
  nn = 0;      % the total number of structures in USPEX_STRUC with current composition
  for i=1:N
    if (USPEX_STRUC.POPULATION(i).numIons == numions)
      nn = nn+1;
      number(nn) = i;
      fitness(nn) = USPEX_STRUC.POPULATION(i).Enthalpies(end);
      enth(nn)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(numions);
    end
  end
  [nothing, ranking] = sort(fitness);

  % Remove all duplicates:
  %--------------------------------------------------------------------------
  GoodList = zeros(max_num, 1); % Good Structure
  item = 1;
  GoodList(1) = 1; % Good Structure
  for i = 2 : nn
    ii = number(ranking(i)); %number in USPEX_STRUC
    if USPEX_STRUC.POPULATION(ii).ToCount > 0
        same = 0;
        for j= 1 : i-1
            jj = number(ranking(j)); %number in USPEX_STRUC
            if USPEX_STRUC.POPULATION(jj).ToCount > 0
                if SameStructure(ii,jj,USPEX_STRUC)
                    USPEX_STRUC.POPULATION(ii).ToCount = 0;
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
  %safesave('../USPEX.mat',USPEX_STRUC); We don't save here

  Good_num = item; %less than max_num
  % output all accepted structures
  string_comp = ['composition_' num2str(comp_x) '_' num2str(comp_y) '_' num2str(comp_z)];
  fpath  = ['__goodStructures/' string_comp];
  fp1 = fopen(fpath, 'w');
  unixCmd(['cat /dev/null > ' fpath]);  %Start to get new POSCAR
  unixCmd(['cat /dev/null > ' fpath '_POSCARS']);
  %fp2 = fopen('goodStructures_POSCARS', 'w');
  fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes   Enthalpies SYMM\n');
  fprintf(fp1,'                       (eV/atom)    (A^3/atom)    (eV)        \n');
  for i = 1 : Good_num
    jj = number(ranking(GoodList(i)));
    num    = USPEX_STRUC.POPULATION(jj).numIons;
    volume = USPEX_STRUC.POPULATION(jj).Vol/sum(num);
    symg   = USPEX_STRUC.POPULATION(jj).symg;
    lattice= USPEX_STRUC.POPULATION(jj).LATTICE;
    coor   = USPEX_STRUC.POPULATION(jj).COORDINATES;
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
        composition=[composition,blanks(shift(length(num)))];
    end
    
    j = ranking(GoodList(i));
    fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %9.4f  %4d\n',...
        jj, composition, enth(j), volume, fitness(j), symg);
    % POSCARS
    Write_POSCAR(atomType, jj, symg, num, lattice, coor);
    unixCmd([' cat POSCAR      >> ' fpath '_POSCARS']);  
  end
  fclose(fp1);

end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------%
%----------------------- case of one-component clusters --------------------------%
%---------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(ORG_STRUC.numIons,2) == 1

for comp = ORG_STRUC.numIons(1,1) : ORG_STRUC.numIons(2,1)
  
  numions = comp;
  fitness = [];
  enth = [];
  number = []; % numner of current structure in USPEX_STRUC (array)
  nn = 0;      % the total number of structures in USPEX_STRUC with current composition
  for i=1:N
    if (USPEX_STRUC.POPULATION(i).numIons == numions)
      nn = nn+1;
      number(nn) = i;
      fitness(nn) = USPEX_STRUC.POPULATION(i).Enthalpies(end);
      enth(nn)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(numions);
    end
  end
  [nothing, ranking] = sort(fitness);

  % Remove all duplicates:
  %--------------------------------------------------------------------------
  GoodList = zeros(max_num, 1); % Good Structure
  item = 1;
  GoodList(1) = 1; % Good Structure
  for i = 2 : nn
    ii = number(ranking(i)); %number in USPEX_STRUC
    if USPEX_STRUC.POPULATION(ii).ToCount > 0
        same = 0;
        for j= 1 : i-1
            jj = number(ranking(j)); %number in USPEX_STRUC
            if USPEX_STRUC.POPULATION(jj).ToCount > 0
                if SameStructure(ii,jj,USPEX_STRUC)
                    USPEX_STRUC.POPULATION(ii).ToCount = 0;
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
  %safesave('../USPEX.mat',USPEX_STRUC); We don't save here

  Good_num = item; %less than max_num
  % output all accepted structures
  string_comp = ['composition_' num2str(comp)];
  fpath  = ['__goodStructures/' string_comp];
  fp1 = fopen(fpath, 'w');
  unixCmd(['cat /dev/null > ' fpath]);  %Start to get new POSCAR
  unixCmd(['cat /dev/null > ' fpath '_POSCARS']);
  %fp2 = fopen('goodStructures_POSCARS', 'w');
  fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes   Enthalpies SYMM\n');
  fprintf(fp1,'                       (eV/atom)    (A^3/atom)    (eV)        \n');
  for i = 1 : Good_num
    jj = number(ranking(GoodList(i)));
    num    = USPEX_STRUC.POPULATION(jj).numIons;
    volume = USPEX_STRUC.POPULATION(jj).Vol/sum(num);
    symg   = USPEX_STRUC.POPULATION(jj).symg;
    lattice= USPEX_STRUC.POPULATION(jj).LATTICE;
    coor   = USPEX_STRUC.POPULATION(jj).COORDINATES;
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
        composition=[composition,blanks(shift(length(num)))];
    end
    
    j = ranking(GoodList(i));
    fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %9.4f  %4d\n',...
        jj, composition, enth(j), volume, fitness(j), symg);
    % POSCARS
    Write_POSCAR(atomType, jj, symg, num, lattice, coor);
    unixCmd([' cat POSCAR      >> ' fpath '_POSCARS']);  
  end
  fclose(fp1);

end

end
