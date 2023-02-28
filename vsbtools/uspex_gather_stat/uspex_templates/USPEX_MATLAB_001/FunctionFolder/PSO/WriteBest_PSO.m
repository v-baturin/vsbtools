function WriteBest_PSO(IND)
global POP_STRUC

fpath = [ POP_STRUC.resFolder '/BESTIndividuals'];
fp = fopen(fpath, 'a+');

gen = POP_STRUC.generation;
Count = POP_STRUC.DoneOrder(IND);
symg = POP_STRUC.POPULATION(IND).symg;
fit = POP_STRUC.POPULATION(IND).FITNESSES(end);
num = POP_STRUC.POPULATION(IND).numIons;

composition = sprintf('%3d',num);
shift=[4, 2, 1]; %so far we only consider 6 component
if size(composition,2)<11
    composition=[composition,blanks(shift(length(num)))];
end

volume = det(POP_STRUC.POPULATION(IND).LATTICE);

if ~isempty(POP_STRUC.POPULATION(IND).K_POINTS)
    KPOINTS =  num2str(POP_STRUC.POPULATION(IND).K_POINTS(end,:));
else
    KPOINTS = ['0 0 0'];
end
if size(KPOINTS,2)<6
    KPOINTS = [KPOINTS,blanks(6-size(KPOINTS,2))];
    KPOINTS = strjust(KPOINTS,'center');
end

entropy = POP_STRUC.struc_entr(IND);
if isempty(entropy)
    entropy = 0;
end

if sum(num)>0
    a_o = sum(POP_STRUC.POPULATION(IND).order)/sum(num);
else
    a_o = 0;
end

if isempty(POP_STRUC.POPULATION(IND).Parents)
    fprintf(fp,'%3d %4d    random    [%11s]  %10.3f   %10.5f  [%7s]  %4s %8.4f %8.4f\n', gen, Count, composition, fit, volume, KPOINTS, symg, entropy, a_o);
elseif isfield(POP_STRUC.POPULATION(IND).Parents, 'mut_degree')
    fprintf(fp,'%3d %4d softmutation [%11s]  %10.3f   %10.5f  [%7s]  %4s %8.4f %8.4f\n', gen, Count, composition, fit, volume, KPOINTS, symg, entropy, a_o);
elseif isfield(POP_STRUC.POPULATION(IND).Parents, 'Seeds')
    fprintf(fp,'%3d %4d    Seeds     [%11s]  %10.3f   %10.5f  [%7s]  %4s %8.4f %8.4f\n', gen, Count, composition, fit, volume, KPOINTS, symg, entropy, a_o);
elseif isfield(POP_STRUC.POPULATION(IND).Parents, 'N_Moved')
    fprintf(fp,'%3d %4d   mutation   [%11s]  %10.3f   %10.5f  [%7s]  %4s %8.4f %8.4f\n', gen, Count, composition, fit, volume, KPOINTS, symg, entropy, a_o);
elseif isfield(POP_STRUC.POPULATION(IND).Parents, 'fracFrac')
    if POP_STRUC.POPULATION(IND).Parents.global==1
        fprintf(fp,'%3d %4d   heredity-g [%11s]  %10.3f   %10.5f  [%7s]  %3s %8.4f %8.4f\n', gen, Count, composition, fit, volume, KPOINTS, symg, entropy, a_o);
    else
        fprintf(fp,'%3d %4d   heredity-l [%11s]  %10.3f   %10.5f  [%7s]  %3s %8.4f %8.4f\n', gen, Count, composition, fit, volume, KPOINTS, symg, entropy, a_o);
    end
elseif isfield(POP_STRUC.POPULATION(IND).Parents, 'keptAsBest')
    fprintf(fp,'%3d %4d   keptBest   [%11s]  %10.3f   %10.5f  [%7s]  %4s %8.4f %8.4f\n', gen, Count, composition, fit, volume, KPOINTS, symg, entropy, a_o);
end

fclose(fp);
