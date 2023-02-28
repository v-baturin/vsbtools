%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PRINTING GRAPHS NEEDED FOR DEBUGGING, FOLDER "graphFolder" %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_graphs_3comp(graphFolder, fitness)
global ORG_STRUC
global POP_STRUC
global CLUSTERS

numions = ORG_STRUC.numIons;
numcomp = CLUSTERS.number_ConvexHall;

minions1 = numions(1,1);
minions2 = numions(1,2);
minions3 = numions(1,3);

nx = numions(2,1) - numions(1,1) + 1;
ny = numions(2,2) - numions(1,2) + 1;
nz = numions(2,3) - numions(1,3) + 1;

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
s = 400/sqrt(nx*ny*nz); %size of circles

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
%%%%%%%%%%%%%%%%%%%%%%%%%% Fitness vs numions %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    h = figure;
    hold;
    set(gcf,'Visible','off'); % Use this switcher to prevent Matlab foregroung printing
    xlabel('numIons');
    ylabel('Fitness');
    
    %%% Line coonnecting fitnesses of best clusters of all compositions %%
    xx = []; ff = [];
    n=0;
    for cx = numions(1,1):numions(2,1)
    for cy = numions(1,2):numions(2,2)
    for cz = numions(1,3):numions(2,3)
        
        cxx = cx-minions1+2;
        cyy = cy-minions2+2;
        czz = cz-minions3+2;
        
        Ex_mean = 0.5*(CLUSTERS.composition(cxx+1,cyy,czz).bestEnthalpy + CLUSTERS.composition(cxx-1,cyy,czz).bestEnthalpy);
        Ey_mean = 0.5*(CLUSTERS.composition(cxx,cyy+1,czz).bestEnthalpy + CLUSTERS.composition(cxx,cyy-1,czz).bestEnthalpy);
        Ez_mean = 0.5*(CLUSTERS.composition(cxx,cyy,czz+1).bestEnthalpy + CLUSTERS.composition(cxx,cyy,czz-1).bestEnthalpy);
        Eref = min([CLUSTERS.composition(cxx,cyy,czz).bestEnthalpy, Ex_mean, Ey_mean, Ez_mean]);

        n=n+1;
        xx(n) = n; %cx + (1/ny) * (cy-numions(1,2));
        ff(n) = CLUSTERS.composition(cxx,cyy,czz).bestEnthalpy - Eref;
    end
    end
    end
    plot(xx,ff,'LineWidth',1);
    
    %%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%%%
    for k = 1:size(operators,1)
        x=[]; F=[]; n=0;
        for i=1:length(POP_STRUC.POPULATION)
            if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, operators(k,:)))
                ix = POP_STRUC.POPULATION(i).numIons(1);
                iy = POP_STRUC.POPULATION(i).numIons(2);
                iz = POP_STRUC.POPULATION(i).numIons(3);
                
                n = n + 1;
                x(n) = (ix-minions1)*ny*nz + (iy-minions2)*nz + (iz-minions3) + 1;
                F(n) = fitness(i);
            end
        end
        scatter(x,F,s,ccc(k,:),'filled','MarkerEdgeColor','k');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% variation operators %%%%%%%%%%%%%%%%%%%%%%%%%%
    min_ax = 1;
    max_ax = nx*ny*nz + 4;
    min_ay = 0;
    max_ay = 10; %CLUSTERS.fitness_range;
%     xxx = max_ax - 0.5/(numions(2,2)-numions(1,2)+1); % - (max_ax-min_ax)/10;
    
    for i = 1:size(operators_text,1)
        yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
        scatter(max_ax-3.5,yyy,60,ccc(i,:),'filled','MarkerEdgeColor','k');
        text(max_ax-3.5+(max_ax-min_ax)/50, yyy, ['(' num2str(num_opera(i)) ') ' operators_text(i,:)], 'FontSize', 10);
    end
    axis([min_ax,max_ax,min_ay,max_ay]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% ticks & labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n=0;
    x_tick = [];
    s_tick = [];
    
    for ix = numions(1,1):numions(2,1)
    for iy = numions(1,2):numions(2,2)
    for iz = numions(1,3):numions(2,3)
        n=n+1;
        x_tick(n) = n;
        s_tick{n} = [num2str(ix) ',' num2str(iy) ',' num2str(iz)];
    end
    end
    end
    
    set(gca, 'xtickLabelRotation', -45);
    set(gca,'XTick', x_tick,'XTickLabel',s_tick);
    print(h,'-dtiff','-r120',[graphFolder '/Fitness_gen' num2str(POP_STRUC.generation) '.tif']);
    close(h);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Energy - referenceEnergy vs numions %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% reading data from the file "Reference_Data.txt" %%%%%%%%%%%%%%

el1 = megaDoof(ORG_STRUC.atomType(1));
el2 = megaDoof(ORG_STRUC.atomType(2));
el3 = megaDoof(ORG_STRUC.atomType(3));

E_ref = [];
filename = ['Reference_Data_' el1 '_' el2 '_' el3 '.txt'];
fid_ref = fopen([ORG_STRUC.homePath '/' filename],'rt');

if (fid_ref~=-1)
    ss = fgets(fid_ref);
    ss = fgets(fid_ref);
    max_at = str2num(ss);
    
    for n = 1:max_at(1)
        ss = fgets(fid_ref);
        if n~= str2num(ss)
            disp('error in reference data!');
        end
        
        for m = 1:max_at(2)
            ss = fgets(fid_ref);
            sn = str2num(ss);
            for k = 1:max_at(3)
                E_ref(n,m,k) = sn(k);
            end
        end
    end
    fclose(fid_ref);
end

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
        n=0;
        for ix = numions(1,1):numions(2,1)
        for iy = numions(1,2):numions(2,2)
        for iz = numions(1,3):numions(2,3)
            ixx = ix-minions1+2;
            iyy = iy-minions2+2;
            izz = iz-minions3+2;
                        
            n = n + 1;
            Xrel(n) = n;
            Erel(n) = CLUSTERS.composition(ixx,iyy,izz).bestEnthalpy - E_ref(ix,iy,iz);
        end
        end
        end
        plot(Xrel,Erel,'LineWidth',1);
        
        %%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%
        Xmag=[]; 
        Emag=[]; 
        n=0;
        for i = 1 : length(POP_STRUC.POPULATION)
            if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, 'convex'))
                ix = POP_STRUC.POPULATION(i).numIons(1);
                iy = POP_STRUC.POPULATION(i).numIons(2);
                iz = POP_STRUC.POPULATION(i).numIons(3);
                
                n = n + 1;
                Xmag(n) = (ix-minions1)*ny*nz + (iy-minions2)*nz + (iz-minions3) + 1;
                Emag(n) = POP_STRUC.POPULATION(i).Enthalpies(end) - E_ref(ix,iy,iz);
            end
        end
        scatter(Xmag,Emag,s,'k','filled','MarkerEdgeColor','k');
        
        %%%%%%%%%%%%%%%%%%%%% diapason of axis %%%%%%%%%%%%%%%%%%%%%%%%%%%
        min_ax = 1;
        max_ax = nx*ny*nz;
        min_ay = -1;
        max_ay = 3;
        axis([min_ax,max_ax,min_ay,max_ay]);
        plot([min_ax,max_ax],[0,0],'k');
        
        %%%%%%%%%%%%%%%%%%%%%%% ticks & labels %%%%%%%%%%%%%%%%%%%%%%%%%%%
        n=0;
        x_tick = [];
        s_tick = [];
        
        for ix = numions(1,1):numions(2,1)
        for iy = numions(1,2):numions(2,2)
        for iz = numions(1,3):numions(2,3)
            n=n+1;
            x_tick(n) = n;
            s_tick{n} = [num2str(ix) ',' num2str(iy) ',' num2str(iz)];
        end
        end
        end
        
        set(gca, 'xtickLabelRotation', -45);
        set(gca,'XTick', x_tick,'XTickLabel',s_tick);
        box on;
        print(h,'-dtiff','-r120',[graphFolder '/Energy-refEnergy_gen' num2str(POP_STRUC.generation) '.tif']);
        close(h);
    catch
    end
end
