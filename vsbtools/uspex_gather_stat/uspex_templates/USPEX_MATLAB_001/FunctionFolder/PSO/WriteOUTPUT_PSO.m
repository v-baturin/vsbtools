function WriteOUTPUT_PSO(Ind_No, resFolder)

%To update all the necessary items after each structure is done
%1, OUTPUT.txt
%2, Individuals
%3, Origin
%Lastly updated by Qiang Zhu (2014/02/18)

global PSO_STRUC
gen      = PSO_STRUC.POPULATION(Ind_No).gen;
symg     = PSO_STRUC.POPULATION(Ind_No).symg;
enth     = PSO_STRUC.POPULATION(Ind_No).Enthalpies(end);
num      = PSO_STRUC.POPULATION(Ind_No).numIons;
fit      = PSO_STRUC.POPULATION(Ind_No).Fitness;
howcome  = PSO_STRUC.POPULATION(Ind_No).howCome;
KPOINTS  = PSO_STRUC.POPULATION(Ind_No).K_POINTS(end,:);
volume   = PSO_STRUC.POPULATION(Ind_No).Vol;
density  = PSO_STRUC.POPULATION(Ind_No).density;
entropy  = PSO_STRUC.POPULATION(Ind_No).struc_entr;
s        = PSO_STRUC.POPULATION(Ind_No).S_order;
order    = PSO_STRUC.POPULATION(Ind_No).order;

composition = sprintf('%3d',num);
shift=[4, 2, 1]; %so far we only consider 6 component
if size(composition,2)<11
    composition=[composition,blanks(shift(length(num)))];
end

if ~isempty(PSO_STRUC.POPULATION(Ind_No).Parents)
    par_ID = PSO_STRUC.POPULATION(Ind_No).Parents.parent;
    par_fit= PSO_STRUC.POPULATION(Ind_No).Parents.enthalpy;
else
    par_ID = '0';
    par_fit = enth/sum(num);
end

if isempty(entropy)
    entropy = 0;
end

if sum(num)>0
    a_o = sum(order)/sum(num);
else
    a_o = 0;
end

fpath =  [resFolder '/OUTPUT.txt'];
fpath1 =  [resFolder '/Individuals'];
fpath2 = [resFolder '/origin'];
fp = fopen(fpath, 'a+');
fp1 = fopen(fpath1, 'a+');
fp2 = fopen(fpath2, 'a+');

fprintf(fp,'%4d %-11s [%11s] %10.3f  %10.3f   [%2d %2d %2d]  %3d\n', Ind_No, howcome, composition, enth, volume, KPOINTS(:), symg);
fprintf(fp1,'%3d %4d %-11s [%11s] %9.3f %9.3f %7.3f     N/A    [%2d %2d %2d] %3d  %6.3f %6.3f %6.3f\n', ...
    gen, Ind_No, howcome, composition, enth, volume, density, KPOINTS(:), symg, entropy, a_o, s);

fprintf(fp2,'%4d %-11s %8.3f  %8.3f  [%10s]\n', Ind_No, howcome, enth/sum(num), par_fit, par_ID);
fclose(fp);
fclose(fp1);
fclose(fp2);

unixCmd(['echo ' num2str(PSO_STRUC.POPULATION(Ind_No).Enthalpies,'%10.3f') ' >> ' resFolder '/enthalpies_complete.dat']);

