function WriteIndividual(resFolder)

%To update all the necessary items at the end of each generation
%Individuals would be updated
%Lastly updated by Qiang Zhu (2014/02/18)

global PSO_STRUC

fpath = [ resFolder '/Individuals'];
fp = fopen(fpath, 'w');
fprintf(fp, 'Gen   ID    Origin   Composition    Enthalpy   Volume  Density   Fitness   KPOINTS  SYMM  Q_entr A_order S_order\n');
fprintf(fp, '                                      (eV)      (A^3)  (g/cm^3)\n');

for i=1:length(PSO_STRUC.POPULATION)
    gen      = PSO_STRUC.POPULATION(i).gen;
    symg     = PSO_STRUC.POPULATION(i).symg;
    enth     = PSO_STRUC.POPULATION(i).Enthalpies(end);
    num      = PSO_STRUC.POPULATION(i).numIons;
    fit      = PSO_STRUC.POPULATION(i).Fitness;
    howcome  = PSO_STRUC.POPULATION(i).howCome;
    KPOINTS  = PSO_STRUC.POPULATION(i).K_POINTS(end,:);
    volume   = PSO_STRUC.POPULATION(i).Vol;
    density  = PSO_STRUC.POPULATION(i).density;
    entropy  = PSO_STRUC.POPULATION(i).struc_entr;
    s        = PSO_STRUC.POPULATION(i).S_order;
    order    = PSO_STRUC.POPULATION(i).order;
    
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
        composition=[composition,blanks(shift(length(num)))];
    end
    
    if isempty(entropy)
        entropy = 0;
    end
    
    if sum(num)>0
        a_o = sum(order)/sum(num);
    else
        a_o = 0;
    end
    fprintf(fp,'%3d %4d %-11s [%11s] %9.3f %9.3f %7.3f %10.3f [%2d %2d %2d] %3d  %6.3f %6.3f %6.3f\n', ...
        gen, i, howcome, composition, enth, volume, density, fit, KPOINTS(:), symg, entropy, a_o, s);
end
fclose(fp);

