%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PRINTING GRAPHS NEEDED FOR DEBUGGING, FOLDER "graphFolder" %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_graphs_1comp(graphFolder, fitness)
global ORG_STRUC
global POP_STRUC
global CLUSTERS

numions = ORG_STRUC.numIons;

%%%%%%%%%% reading data from the file "Reference_Data.txt" %%%%%%%%%%%%%%%
fid_ref = fopen([ORG_STRUC.homePath '/Reference_Data.txt'],'rt');
if (fid_ref~=-1)
    Nmin = fscanf(fid_ref,'%f',1);
    Nmax = fscanf(fid_ref,'%f',1);
    ss = fgets(fid_ref);
    for i = Nmin:Nmax
        N_atoms = fscanf(fid_ref,'%f',1);
        E_RefData(N_atoms) = fscanf(fid_ref,'%f',1);
        ss = fgets(fid_ref);
    end
    fclose(fid_ref);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% colors of operators %%%%%%%%%%%%%%%%%%%%%%%%%%
operators = ['Random'; 'Heredi'; 'H1ered'; 'softmu'; ...
             'CoorMu'; 'Permut'; 'TransM'; ...
             'AddAto'; 'Remove'; ...
             'keptBe'; 'convex'; ' Seeds'];
ccc = [1 1 1; 1 0 0; 1 0.6 0.6; 0 1 0;
       0.3 0.7 0.4; 0.9 0.6 0.9; 1 0 1;
       1 1 0; 0.1 0.5 1;
       0 1 1; 0 0 0; 0.7 0.7 0.7];
operators_text = ['Random       '; 'Heredity     '; 'H1eredity    '; 'Softmutation '; ...
                  'Coormutation '; 'Permutation  '; 'Transmutation'; ...
                  'AddAtom      '; 'RemoveAtom   '; ...
                  'keptBest     '; 'ConvexHull   '; 'Seeds        '];
s = 1500/CLUSTERS.number_compositions; %size of circles

%%%%%%%%%%%% number of structures produced by each operator %%%%%%%%%%%%%%
num_opera = [];
for k = 1 : size(operators,1)
    num_opera(k) = 0;
    for i = 1 : length(POP_STRUC.POPULATION)
        if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, operators(k,:)))
            num_opera(k) = num_opera(k) + 1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Relative energies vs numIons %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    h = figure;
    hold;
    set(gcf,'Visible','off');
    xlabel('numIons');
    ylabel('Relative energy');
    
    N1 = numions(1,1);
    N2 = numions(2,1);
    
    %%%%%%%%%%%%%%%%%%% Plotting reference data %%%%%%%%%%%%%%%%%%%%%%%%%%
    if (fid_ref~=-1)
        x_ref = []; 
        E_ref = [];
        for i = N1:N2
            E_ref(i) = E_RefData(i) - (E_RefData(N1)+(i-N1)*(E_RefData(N2)-E_RefData(N1))/(N2-N1));
            plot([i-1/2, i+1/2], [E_ref(i),E_ref(i)],'k','LineWidth',1);
        end
        for i = 1:N1-1
            E_ref(i)=E_ref(N1);
        end
    end
    
    %%% plotting line connecting the best clusters of all compositions %%%
    x_best = []; 
    E_best = [];
    if (fid_ref~=-1)
        E_Ref1 = E_RefData(N1);
        E_Ref2 = E_RefData(N2);
    else
        E_Ref1 = CLUSTERS.composition(1).bestEnthalpy;
        E_Ref2 = CLUSTERS.composition(N2-N1+1).bestEnthalpy;
    end
    for i = N1:N2
        x_best(i) = i;
        E_best(i) = CLUSTERS.composition(i-N1+1).bestEnthalpy - (E_Ref1+(i-N1)*(E_Ref2-E_Ref1)/(N2-N1));
    end
    for i = 1:N1-1
        E_best(i)=E_best(N1);
    end
    plot(x_best,E_best,'--','LineWidth',1);
    
    %%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%%%
    for k = 1 : size(operators,1)
        x=[]; 
        E=[]; 
        n=0;
        for i=1 : length(POP_STRUC.POPULATION)
            if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, operators(k,:)))
                n = n + 1;
                x(n) = POP_STRUC.POPULATION(i).numIons;
                E(n) = POP_STRUC.POPULATION(i).Enthalpies(end) - (E_Ref1+(x(n)-N1)*(E_Ref2-E_Ref1)/(N2-N1));
            end
        end
        scatter(x,E,s,ccc(k,:),'filled','MarkerEdgeColor','k');
    end
    
    %%%%%%%%%%%%%%% variation operators, diapason of axis %%%%%%%%%%%%%%%%
    E = [];
    for i = 1 : length(POP_STRUC.POPULATION)
        E(i) = POP_STRUC.POPULATION(i).Enthalpies(end) - (E_Ref1+(POP_STRUC.POPULATION(i).numIons-N1)*(E_Ref2-E_Ref1)/(N2-N1));
    end
    min_ax = numions(1,1);
    max_ax = numions(2,1);
    min_ay = min(E_best) - (max(E_best)-min(E_best))/10;
    max_ay = max(E_best) + 3*(max(E_best)-min(E_best));
    %min_ay = min(E_best)-std(E)/20;
    %max_ay = min(E_best)+2*std(E);
    if (fid_ref~=-1)
        %  min_ay = min(min(E_ref),min(E_best))-std(E)/20;
        %  max_ay = min(min(E_ref),min(E_best))+2*std(E);
        min_ay = min(min(E_ref),min(E_best))-(max(E_best)-min(E_best))/10;
        max_ay = max(max(E_ref),max(E_best))+3*(max(E_best)-min(E_best));
    end
    
    xxx = max_ax - (max_ax-min_ax)/10;
    for i = 1:size(operators_text,1);
        yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
        scatter(xxx,yyy,60,ccc(i,:),'filled','MarkerEdgeColor','k');
        text(xxx+(max_ax-min_ax)/50, yyy, ['(' num2str(num_opera(i)) ') ' operators_text(i,:)]);
    end
    if (fid_ref~=-1)
        plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
        text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
    end
    
    axis([min_ax,max_ax,min_ay,max_ay]);
    set(gca,'XTick', N1:N2)
    print(h,'-dtiff','-r120',[graphFolder '/Relative_energies_gen' num2str(POP_STRUC.generation) '.tif']);
    close(h);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Fittnes vs numIons %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    h = figure;
    hold;
    set(gcf,'Visible','off');
    xlabel('numIons');
    ylabel('Fitness');
    for k = 1 : size(operators,1);
        x = []; 
        F = []; 
        n = 0;
        for i = 1 : length(POP_STRUC.POPULATION)
            if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, operators(k,:)))
                n = n + 1;
                x(n) = POP_STRUC.POPULATION(i).numIons;
                F(n) = fitness(i);
            end
        end
        scatter(x,F,s,ccc(k,:),'filled','MarkerEdgeColor','k');
    end
    
    %%%%%%%%% Fitnesses of best clusters in current generation %%%%%%%%%%%
    xx = []; 
    FF = [];
    for i = N1:N2
        xx(i) = i;
        ii = 1;
        while 1
            ii = ii + 1;
            if (CLUSTERS.ConvexHall(ii).numIons > i)
                break;
            end
        end
        nn1 = CLUSTERS.ConvexHall(ii-1).numIons;
        nn2 = CLUSTERS.ConvexHall(ii).numIons;
        EE1 = CLUSTERS.ConvexHall(ii-1).Enthalpy;
        EE2 = CLUSTERS.ConvexHall(ii).Enthalpy;
        FF(i) = (CLUSTERS.composition(i-N1+1).bestEnthalpy - (EE1+(i-nn1)*(EE2-EE1)/(nn2-nn1))); % / i;
    end
    plot(xx,FF,'--','LineWidth',1);
    
    %%%%%%%%%%%%%%%% variation operators, diapason of axis %%%%%%%%%%%%%%%
    min_ax = numions(1,1);
    max_ax = numions(2,1);
    min_ay = 0;
    max_ay = 5; %CLUSTERS.fitness_range;
    xxx = max_ax - (max_ax-min_ax)/10;
    for i = 1 : size(operators_text,1);
        yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
        scatter(xxx,yyy,60,ccc(i,:),'filled','MarkerEdgeColor','k');
        text(xxx+(max_ax-min_ax)/50, yyy, ['(' num2str(num_opera(i)) ') ' operators_text(i,:)]);
    end
    if (fid_ref~=-1)
        plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
        text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
    end
    
    axis([min_ax,max_ax,min_ay,max_ay]);
    set(gca,'XTick', N1:N2)
    print(h,'-dtiff','-r120',[graphFolder '/Fitness_gen' num2str(POP_STRUC.generation) '.tif']);
    close(h);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Energy/numions vs numions %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    h = figure;
    hold;
    set(gcf,'Visible','off');
    xlabel('numIons');
    ylabel('Energy / numIons');
    
    %%%%%%%%%%%%%%%%%%%%%%% Plotting reference data %%%%%%%%%%%%%%%%%%%%%%
    if (fid_ref~=-1)
        E_ref = [];
        for i = 1:N1-1
            E_ref(i) = 10000;
        end
        for i = N1:N2
            E_ref(i) = E_RefData(i) / i;
            plot([i-1/2, i+1/2], [E_ref(i),E_ref(i)],'k','LineWidth',1);
        end
    end
    
    %%% plotting line connecting the best clusters of all compositions %%%
    x_best = [];
    E_best = [];
    for i = N1:N2
        x_best(i) = i;
        E_best(i) = CLUSTERS.composition(i-N1+1).bestEnthalpy / i;
    end
    for i = 1 : N1-1
        E_best(i) = E_best(N1);
    end
    plot(x_best,E_best,'--','LineWidth',1);
    
    %%%%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%
    for k = 1 : size(operators,1)
        x = []; 
        E = []; 
        n = 0;
        for i = 1 : length(POP_STRUC.POPULATION)
            if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, operators(k,:)))
                n = n + 1;
                x(n) = POP_STRUC.POPULATION(i).numIons;
                E(n) = POP_STRUC.POPULATION(i).Enthalpies(end) / x(n);
            end
        end
        scatter(x,E,s,ccc(k,:),'filled','MarkerEdgeColor','k');
    end
    
    %%%%%%%%%%%%%%%%% variation operators, diapason of axis %%%%%%%%%%%%%%
    E = [];
    for i = 1 : length(POP_STRUC.POPULATION)
        E(i) = POP_STRUC.POPULATION(i).Enthalpies(end) / POP_STRUC.POPULATION(i).numIons;
    end
    min_ax = numions(1,1);
    max_ax = numions(2,1);
    min_ay = min(E_best) - (max(E_best)-min(E_best))/20;
    max_ay = max(E_best) + (max(E_best)-min(E_best))/2;
    if (fid_ref~=-1)
        min_ay = min(min(E_ref),min(E_best))-(max(E_best)-min(E_best))/10;
    end
    
    xxx = max_ax - (max_ax-min_ax)/10;
    for i = 1:size(operators_text,1);
        yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
        scatter(xxx,yyy,60,ccc(i,:),'filled','MarkerEdgeColor','k');
        text(xxx+(max_ax-min_ax)/50, yyy, ['(' num2str(num_opera(i)) ') ' operators_text(i,:)]);
    end
    if (fid_ref~=-1)
        plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
        text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
    end
    
    axis([min_ax,max_ax,min_ay,max_ay]);
    set(gca,'XTick', N1:N2)
    print(h,'-dtiff','-r120',[graphFolder '/Energy_div_numIons_gen' num2str(POP_STRUC.generation) '.tif']);
    close(h);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Energy - referenceEnergy vs numions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (fid_ref~=-1)
    try
        h = figure;
        hold;
        set(gcf,'Visible','off'); % Use this switcher to prevent Matlab foregroung printing
        xlabel('numIons');
        ylabel('Energy - referenceEnergy');
        
        % Plotting line connecting best clusters of all compositions in current generations %
        Erel = [];
        Xrel = [];
        n = 0;
        for i = numions(1,1):numions(2,1)
            n = n + 1;
            Erel(n) = CLUSTERS.composition(i-numions(1,1)+1).bestEnthalpy - E_RefData(i);
            Xrel(n) = i;
        end
        plot(Xrel,Erel,'LineWidth',1);
        dlmwrite([graphFolder '/Energy-refData_gen' num2str(POP_STRUC.generation) '.txt'],Erel');

        %%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%
        Xmag = []; 
        Emag = []; 
        n = 0;
        for i = 1 : length(POP_STRUC.POPULATION)
            if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, 'convex'))
                n = n + 1;
                Xmag(n) = POP_STRUC.POPULATION(i).numIons;
                Emag(n) = POP_STRUC.POPULATION(i).Enthalpies(end) - E_RefData(Xmag(n));
            end
        end
        scatter(Xmag,Emag,s,'k','filled','MarkerEdgeColor','k');
        
        %%%%%%%%%%%%%%%%%%%%%%% diapason of axis %%%%%%%%%%%%%%%%%%%%%%%%%
        min_ax = numions(1,1);
        max_ax = numions(2,1);
        min_ay = -0.5;
        max_ay = 3;
        axis([min_ax,max_ax,min_ay,max_ay]);
        plot([min_ax,max_ax],[0,0],'k');
        
        print(h,'-dtiff','-r120',[graphFolder '/Energy-refEnergy_gen' num2str(POP_STRUC.generation) '.tif']);
        close(h);
    catch
    end
end
