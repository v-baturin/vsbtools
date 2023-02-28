function makeFigures(pickUpNCount, Nsteps, ConstLat)

%All the data are loaded from USPEX_STRUC
%Lastly updated by Qiang Zhu (2014/02/18)

global USPEX_STRUC

if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end

N         = length(USPEX_STRUC.POPULATION);
N_GEN     = length(USPEX_STRUC.GENERATION);
enth      = zeros(N,Nsteps);
vol       = zeros(N,1);
fit       = zeros(N,1);
ID        = zeros(N,1);

for i=1:N
    numIons   = USPEX_STRUC.POPULATION(i).numIons;
    fit(i)    = USPEX_STRUC.POPULATION(i).Fitness;
    enth(i,:) = USPEX_STRUC.POPULATION(i).Enthalpies/sum(numIons);
    vol(i)    = det(USPEX_STRUC.POPULATION(i).LATTICE)/sum(numIons);
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
    delete('Energy_vs_N.pdf');
    x = []; y = [];
    x = ID;
    y = enth(:,Nsteps);
    h = figure;
    set(gcf,'Visible','off');% Use this switcher to prevent Matlab foreground 
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    delete('Fitness_vs_N.pdf');
    x = []; y = [];
    x = ID;
    y = fit;
    y_max = max(y(find(y<1000))); %kick out all the structures with fitness=10000;
    y(find(y>=10000))=y_max;
    h = figure;
    set(gcf,'Visible','off');    % Use this switcher to prevent Matlab foreground
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ConstLat~=1
    try
        delete('Energy_vs_Volume.pdf');
        x = []; y = [];
        x = enth(:,Nsteps);
        y = vol;
        h = figure;
        set(gcf,'Visible','off'); % Use this switcher to prevent Matlab foreground
        if OctaveMode == 0
            scatter(x,y,'MarkerEdgeColor','k','MarkerFaceColor','g');
        else
            scatter(x,y);
        end
        xlabel('Enthalpy');
        ylabel('Volume');
        print(h,'-dpdf','Energy_vs_Volume.pdf');
    catch
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphs En(En-1) - energy change between different optimisation steps
try
    e_complete = enth;
    sub = ceil((Nsteps-1)/2);
    delete('E_series.pdf');
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

try
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1 (Echild - Eparent(s)) as function of N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delete('Variation-Operators.pdf');
    h = figure;
    set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foreground
    % 0 - latmutation,  1 - softmutation & coormutation, 2 - rotation,
    % 3 - permutation,  4 - heredity
    N_count = N-pickUpNCount;
    x = ID(N-N_count+1:N);
    y = zeros(N_count,1);
    c = zeros(N_count,1);
    s = ones(N_count,1);
    handle = fopen('origin');
    tmp = fgetl(handle);
    for i=1:N_count
        tmp = fgetl(handle);
        if tmp ~= -1
            if ~isempty(findstr(tmp, 'LatMutated'))
                c(i) = 0;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'CoorMutate')) %merged in softmutation!!!
                c(i) = 1;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'softmutate'))
                c(i) = 1;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Rotated'))
                c(i) = 2;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Permutate'))
                c(i) = 3;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Heredity'))
                c(i) = 4;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'Spinmutate'))
                c(i) = 5;
                s(i) = 30;
            elseif ~isempty(findstr(tmp, 'TransMutate'))
                c(i) = 6;
                s(i) = 30;
            end
            CE = str2num(tmp(17:26));
            PE = str2num(tmp(27:36));
            y(i) = CE-PE;
        end
    end
    status = fclose(handle);
    subplot(2,1,1);
    if OctaveMode == 0
        scatter(x,y,s,c,'filled','MarkerEdgeColor','k');
    else
        scatter(x,y,s,c);
    end
    xlabel('Structure index');
    ylabel('E_{Child} - E_{parent}');
    hcb = colorbar('YTick',[0:1:6],'YTickLabel',...
    {'Latmut','Atommut','Rotmut','Permut','Heredity','Spinmut','Transmut'});
    set(hcb,'YTickMode','manual')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fraction versus N
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[a, b] = unix(['grep frac OUTPUT.txt |cut -c56-60']);
    %numOperater = 8;
    if N_GEN > 2
        Parent =  {'Random',        'RandTop',  'Heredity', 'H1eredity', ...
                   'Permutate', 'TransMutate', 'LatMutate',...
                   'softmutate',     'Rotate',      'Spin',...
                   'SecSwitch', 'ShiftBorder','AddAtom','RemoveAtom'};
        
        fractions = {'fracRand', 'fracRandTop',  'fracGene', 'fracG1ene', ...
                     'fracPerm', 'fracTrans',  'fracLatMut',...
                     'fracAtomsMut','fracRotMut','fracSpin',...
                     'fracSecSwitch', 'fracShiftBorder','fracAddAtom','fracRemAtom'};
        color = {'g', 'k', 'k', 'k--', '-g', 'b', 'k--', ...
                 'c-*', 'y-', 'b--', 'y-', 'b--','k-^','k-v'};
        numOperation = length(Parent);
        %str = {'Heredity', 'Random', 'TR Random', 'Softmut',...
        %'Permut', 'Latmut', 'Rotmut', 'Spinmut', 'Transmut'};
        subplot(2,1,2);
        hold on;
        str = '';
        for i = 1:numOperation
            Y = zeros(N_GEN,1);
            for j = 1:N_GEN
                %disp([ 'Y(j) = USPEX_STRUC.GENERATION(j).' fractions{i} ';' ]);
                eval([ 'Y(j) = USPEX_STRUC.GENERATION(j).' fractions{i} ';' ]);
            end

            if sum(Y) > 0
               plot(Y,color{i});
               str{end+1} = Parent{i};
            end
        end

        legend(str);
        xlabel('Generation number');
        ylabel('Fraction of each variation');
        axis([1 N_GEN 0 1]);
    end
    
    print(h,'-dpdf', 'Variation-Operators.pdf');
catch
end
