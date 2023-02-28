function makeFigures(pickUpGen, Nsteps, ConstLat)
%All the data are loaded from PSO_STRUC
%Lastly updated by Qiang Zhu (2014/02/18)

global PSO_STRUC

if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end

N         = length(PSO_STRUC.POPULATION);
enth      = zeros(N,Nsteps);
vol       = zeros(N,1);
fit       = zeros(N,1);
ID        = zeros(N,1);

for i=1:N
    numIons   = PSO_STRUC.POPULATION(i).numIons;
    fit(i)    = PSO_STRUC.POPULATION(i).Fitness;
    enth(i,:) = PSO_STRUC.POPULATION(i).Enthalpies/sum(numIons);
    vol(i) =det(PSO_STRUC.POPULATION(i).LATTICE)/sum(numIons);
    ID(i) = i;
end


% creates a set of figures in pdf format
% Nsteps = number of optimization steps
% E(n), fitness(n), E(V), Echild(Eparents), En(En-1) where n is optimization step
% if ConstLat = 1, we don't output E(V).
% Octave is not able to plot this figure, so let's skip it
% Last updated by Qiang Zhu (2014/02/20)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    unixCmd(['rm Energy_vs_N.pdf']);
    x = []; y = [];
    x = ID;
    y = enth(:,Nsteps);
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    if OctaveMode == 0
        scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
    else
        scatter(x,y);
    end
    xlabel('Structure number');
    ylabel('Enthalpy');
    print(h,'-dpdf', 'Energy_vs_N.pdf');
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    unixCmd(['rm Fitness_vs_N.pdf']);
    x = []; y = [];
    x = ID;
    y = fit;
    y_max = max(y(find(y<1000))); % this is to kick out all the structures with fitness=10000;
    y(find(y>=10000))=y_max;
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
    if OctaveMode == 0
        scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
    else
        scatter(x,y);
    end
    xlabel('Structure number');
    ylabel('Fitness');
    print(h,'-dpdf','-r120','Fitness_vs_N.pdf');
catch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ConstLat~=1
    try
        unixCmd(['rm Energy_vs_Volume.pdf']);
        x = []; y = [];
        x = enth(:,Nsteps);
        y = vol;
        h = figure;
        set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
        if OctaveMode == 0
            scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
        else
            scatter(x,y);
        end
        xlabel('Volume');
        ylabel('Enthalpy');
        print(h,'-dpdf','Energy_vs_Volume.pdf');
    catch
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs En(En-1) - to compare the energy change between different optimisation steps
try
    e_complete = enth;
    sub = ceil((Nsteps-1)/2);
    unixCmd(['rm E_series.pdf']);
    for g = 1 : Nsteps - 1
        x = []; y = [];
        x = e_complete(:,g);
        y = e_complete(:,g+1);
        subplot(sub,2,g);
        if OctaveMode == 0
            scatter(x,y,10,'MarkerEdgeColor','k','MarkerFaceColor','g');
        else
            scatter(x,y,10);
        end
        xlabel(['E' num2str(g)]);
        ylabel(['E' num2str(g+1)]);
    end
    print(h,'-dpdf', 'E_series.pdf');
catch
end
