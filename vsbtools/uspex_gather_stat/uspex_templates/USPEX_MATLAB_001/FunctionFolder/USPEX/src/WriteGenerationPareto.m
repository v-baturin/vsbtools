function WriteGenerationPareto(Ranking, ParetoFront, fitness, Enthalpies, Properties, paretoRanking, opt_sign, ERROR)
global POP_STRUC
resFolder  = POP_STRUC.resFolder;
if paretoRanking > 1
    Properties = -1 * Properties;
end
fitness = opt_sign * -1 * fitness;

fpath = [ resFolder '/ParetoRankingInGeneration'];
fp    = fopen(fpath, 'a+');
fprintf(fp, '---- generation%3d ----\n', POP_STRUC.generation);
if ERROR == 1
    fprintf(fp, 'Good structures cannot satisfy hardness > 2GPa condition, ranking is performed without this constrain! \n');
end
if paretoRanking == 1
    fprintf(fp, 'ParetoFront    ID      Fitness        Convex_hull  \n');
else
    fprintf(fp, 'ParetoFront    ID      Fitness        Convex_hull    3th_Property \n');
end
for i = 1 : length(fitness)
    if paretoRanking == 1 
        fprintf(fp, ' %5d   %8d    %9.3f    %12.3f  \n', ParetoFront(i), Ranking(i), fitness(Ranking(i)), Enthalpies(Ranking(i)));
    else
	fprintf(fp, ' %5d   %8d    %9.3f    %12.3f    %12.3f \n', ParetoFront(i), Ranking(i), fitness(Ranking(i)), Enthalpies(Ranking(i)), Properties(Ranking(i)));
    end
end
fclose(fp);
