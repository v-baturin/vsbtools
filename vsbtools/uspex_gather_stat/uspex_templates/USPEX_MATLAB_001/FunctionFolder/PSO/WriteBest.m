function WriteBest(resFolder)
%Output the BESTINDIVIDUALS info
%Lastly updated by Qiang Zhu (2014/02/18)
global PSO_STRUC
fpath = [ resFolder '/BESTIndividuals'];
fp = fopen(fpath, 'w');
fprintf(fp, 'Gen   ID    Origin   Composition    Enthalpy   Volume  Density   Fitness   KPOINTS  SYMM  Q_entr A_order S_order\n');
fprintf(fp, '                                      (eV)      (A^3)  (g/cm^3)\n');
for i = 1:  length(PSO_STRUC.GENERATION)
    for j=1:length(PSO_STRUC.GENERATION(i).BestID)
        IND     = PSO_STRUC.GENERATION(i).BestID(j);
        gen     = PSO_STRUC.POPULATION(IND).gen;
        fit     = PSO_STRUC.POPULATION(IND).Fitness;
        symg    = PSO_STRUC.POPULATION(IND).symg;
        enth    = PSO_STRUC.POPULATION(IND).Enthalpies(end);
        num     = PSO_STRUC.POPULATION(IND).numIons;
        entropy = PSO_STRUC.POPULATION(IND).struc_entr;
        KPOINTS = PSO_STRUC.POPULATION(IND).K_POINTS(end,:);
        howcome = PSO_STRUC.POPULATION(IND).howCome;
        order   = PSO_STRUC.POPULATION(IND).order;
        s       = PSO_STRUC.POPULATION(IND).S_order;
        volume  = PSO_STRUC.POPULATION(IND).Vol;
        density = PSO_STRUC.POPULATION(IND).density;
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
            gen, IND, howcome, composition, enth, volume, density, fit, KPOINTS(:), symg, entropy, a_o, s);
    end
end
fclose(fp);
