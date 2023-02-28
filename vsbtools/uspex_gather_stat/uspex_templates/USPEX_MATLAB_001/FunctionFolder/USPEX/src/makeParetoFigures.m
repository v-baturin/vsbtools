function makeParetoFigures(fitness, enthalpy, resfolder)

if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end
try
    delete([resfolder '/Fitness_vs_Energy.pdf']);
    x = []; y = [];
    x = enthalpy;
    y = fitness;
    h = figure;
    set(gcf,'Visible','off');% Use this switcher to prevent Matlab foreground 
    if OctaveMode == 0
        scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
    else
        scatter(x,y);
    end
    ylabel('Fitness');
    xlabel('Enthalpy');
    print(h,'-dpdf' , [resfolder '/Fitness_vs_Energy.pdf']);
catch
end
