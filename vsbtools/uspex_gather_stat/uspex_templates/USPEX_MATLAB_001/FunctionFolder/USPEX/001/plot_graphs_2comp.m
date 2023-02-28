%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PRINTING GRAPHS NEEDED FOR DEBUGGING, FOLDER "graphFolder" %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_graphs_2comp(graphFolder, fitness)
global ORG_STRUC
global POP_STRUC
global CLUSTERS

numions = ORG_STRUC.numIons;
numcomp = CLUSTERS.number_ConvexHall;
nx = numions(2,1) - numions(1,1) + 1;
ny = numions(2,2) - numions(1,2) + 1;

%%%%%%%%%%% reading data from the file "Reference_Data.txt" %%%%%%%%%%%%%%
fid_ref = fopen([ORG_STRUC.homePath '/Reference_Data.txt'],'rt');

if (fid_ref~=-1)
    N1_ref = fscanf(fid_ref,'%f',1);
    N2_ref = fscanf(fid_ref,'%f',1);
    ss = fgets(fid_ref);
    
    for i=1:N1_ref
        for j=1:N2_ref
            E_LJ(i,j) = fscanf(fid_ref,'%f',1);
        end
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
s = 400/sqrt(nx*ny); %size of circles

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
%%%%%%%%%%% plotting 3d graph of energies and reference data %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try
%     h = figure;
%     hold;
%     set(gcf,'Visible','off'); % Use this switcher to prevent Matlab foregroung printing
%     ylabel('Energy');
%     
%     %%%%%%%%%% Calculation relative energy of each composition %%%%%%%%%%%
%     x1 = numions(1,1);
%     y1 = numions(1,2);
%     z1 = CLUSTERS.composition(2,2).bestEnthalpy;
%     x2 = numions(2,1);
%     y2 = numions(1,2);
%     z2 = CLUSTERS.composition(x2-x1+2,y2-y1+2).bestEnthalpy;
%     x3 = numions(1,1);
%     y3 = numions(2,2);
%     z3 = CLUSTERS.composition(2,y3-y1+2).bestEnthalpy;
%     
%     E_rel=[];
%     E_rel1=[];
%     n=0;
%     for i = numions(1,1) : numions(2,1)
%         for j = numions(1,2) : numions(2,2)
%             ii = i-numions(1,1)+2;
%             jj = j-numions(1,2)+2;
%             E_rel(ii,jj)=(det([x1 y1 z1;x2 y2 z2;x3 y3 z3]) - ii*det([y1 z1 1;y2 z2 1;y3 z3 1]) + jj*det([x1 z1 1;x2 z2 1;x3 z3 1]))/det([x1 y1 1;x2 y2 1;x3 y3 1]);
%             n=n+1;
%             E_rel1(n)=E_rel(ii,jj);
%         end
%     end
%     
%     delta1=0;
%     for i = numions(1,1) : numions(2,1)
%         for j = numions(1,2) : numions(2,2)
%             ii = i-numions(1,1)+2;
%             jj = j-numions(1,2)+2;
%             delta1 = delta1 + abs(CLUSTERS.composition(ii,jj).bestEnthalpy - E_rel(ii,jj));
%         end
%     end
%     
%     deltaE = delta1/(nx*ny*(ny-1));
%     %deltaE=0;
%     
%     %%%%%%% plotting grid of best structures at current generation %%%%%%%
%     for i = numions(1,1) : numions(2,1)
%         xbest=[]; ybest=[];
%         nn=0;
%         for j= numions(1,2) : numions(2,2)
%             ii = i-numions(1,1)+2;
%             jj = j-numions(1,2)+2;
%             nn = nn + 1;
%             xbest(nn) = ii + numions(1,1) - 2 + (1/ny) * (jj-2);
%             ybest(nn) = CLUSTERS.composition(ii,jj).bestEnthalpy  - E_rel(ii,jj) - (jj-2)*deltaE;
%         end;
%         plot(xbest,ybest,'LineWidth',1);
%     end;
%     
%     for j= numions(1,2) : numions(2,2)
%         xbest=[]; ybest=[];
%         nn=0;
%         for i = numions(1,1) : numions(2,1)
%             ii = i-numions(1,1)+2;
%             jj = j-numions(1,2)+2;
%             nn = nn + 1;
%             xbest(nn) = ii + numions(1,1) - 2 + (1/ny) * (jj-2);
%             ybest(nn) = CLUSTERS.composition(ii,jj).bestEnthalpy  - E_rel(ii,jj) - (jj-2)*deltaE;
%         end;
%         plot(xbest,ybest,'LineWidth',1);
%     end;
%     
%     %%%%%%%%%%%%%%%%%% plotting axes, ticks and labels %%%%%%%%%%%%%%%%%%%
%     if 1==0
%         set(gca,'XTick', 0);
%         set(gca,'YTick', round(min(E_rel1)+1):round(max(E_rel1)+3));
%         text(numions(2,1)+0.5,min(E_rel1)-0.2,'numions (A)');
%         plot([numions(1,1), numions(1,1)+1], [min(E_rel1), min(E_rel1) - (ny+1)*deltaE], 'k', 'LineWidth', 1);
%         plot([numions(1,1), numions(2,1)+1], [min(E_rel1), min(E_rel1)], 'k', 'LineWidth', 1);
%         for i = numions(1,2)+1: numions(2,2)
%             ttx = numions(1,1)+(i-numions(1,2))/ny;
%             tty = min(E_rel1)-(i-numions(1,2))*deltaE*(ny+1)/ny;
%             text(ttx-0.02,tty - 0.15, num2str(i));
%             plot([ttx,ttx] , [tty,tty+0.1],'k', 'LineWidth', 1);
%         end
%         for i = numions(1,1)+1: numions(2,1)
%             text(i-0.02,min(E_rel1)-0.15, num2str(i));
%             plot([i,i], [min(E_rel1),min(E_rel1)+0.1], 'k', 'LineWidth', 1);
%         end
%         text(numions(1,1),min(E_rel1)-0.2, [num2str(numions(1,1)) ',' num2str(numions(1,2))]);
%         text(numions(1,1)+1, min(E_rel1) - (ny+1)*deltaE - 0.1, 'numions (B)');
%         text(numions(1,1)+nx/2-1, max(E_rel1)+2.6, ['AxBy      GENERATION ' num2str(POP_STRUC.generation)]);
%     end
%     
%     %%%%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%
%     for k = 1:size(operators,1)
%         x=[]; E=[]; n=0;
%         for i=1:length(POP_STRUC.POPULATION)
%             if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, operators(k,:)))
%                 ix = POP_STRUC.POPULATION(i).numIons(1);
%                 iy = POP_STRUC.POPULATION(i).numIons(2);
%                 n = n + 1;
%                 x(n) = ix + (1/ny) * (iy-numions(1,2));
%                 E(n) = POP_STRUC.POPULATION(i).Enthalpies(end) - E_rel(ix-numions(1,1)+2,iy-numions(1,2)+2) - (iy-numions(1,2))*deltaE;
%             end
%         end
%         scatter(x,E,s,ccc(k,:),'filled','MarkerEdgeColor','k');
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%% plotting reference data %%%%%%%%%%%%%%%%%%%%%%
%     if (fid_ref~=-1)
%         xref=[]; Eref=[]; n=0;
%         for i = numions(1,1) : numions(2,1)
%             for j = numions(1,2) : numions(2,2)
%                 if E_LJ(i,j)~=0
%                     n=n+1;
%                     xref(n) = i + (1/ny) * (j-numions(1,2));
%                     Eref(n) = E_LJ(i,j)  - E_rel(i-numions(1,1)+2,j-numions(1,2)+2) - (j-numions(1,2))*deltaE;
%                 end
%             end
%         end
%         scatter(xref,Eref,s/3,'k','filled','MarkerEdgeColor','k');
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%% diapason of axis %%%%%%%%%%%%%%%%%%%%%%%%%%
%     min_ax = numions(1,1);
%     max_ax = numions(2,1)+1;
%     nn=0; Ebest = [];
%     for i = numions(1,1) : numions(2,1)
%         for j = numions(1,2) : numions(2,2)
%             ii = i-numions(1,1)+2;
%             jj = j-numions(1,2)+2;
%             nn = nn + 1;
%             Ebest(nn) = CLUSTERS.composition(ii,jj).bestEnthalpy  - E_rel(ii,jj) - (jj-2)*deltaE;
%         end
%     end
%     min_ay = min(Ebest)-(max(Ebest)-min(Ebest))/30;
%     max_ay = max(Ebest) + (max(Ebest)-min(Ebest))/2;
%     if (fid_ref~=-1) && (min(Eref)<min(Ebest))
%         min_ay = min(Eref)-(max(Ebest)-min(Eref))/30;
%         max_ay = max(Ebest) + (max(Ebest)-min(Eref))/2;
%     end
%     axis([min_ax,max_ax,min_ay,max_ay]);
%     
%     %%%%%%%%%%%%%%%%%%%%%%% variation operators %%%%%%%%%%%%%%%%%%%%%%%%%%
%     xxx = max_ax - 0.5/(numions(2,2)-numions(1,2)+1); % - (max_ax-min_ax)/10;
%     for i = 1:size(operators_text,1)
%         yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
%         scatter(xxx,yyy,60,ccc(i,:),'filled','MarkerEdgeColor','k');
%         text(xxx+(max_ax-min_ax)/50, yyy, ['(' num2str(num_opera(i)) ') ' operators_text(i,:)]);
%     end
%     if (fid_ref~=-1)
%         scatter(xxx, yyy-(max_ay-min_ay)/30,60/3,'k','filled','MarkerEdgeColor','k');
%         text(xxx+(max_ax-min_ax)/50, yyy-(max_ay-min_ay)/30, 'Ref. data');
%     end
%     
%     print(h,'-dtiff','-r120',[graphFolder '/Energies_3d_gen' num2str(POP_STRUC.generation) '.tif']);
%     close(h);
% catch
% end

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
            done = 0;
            for i = 1:numcomp
                x_ch = CLUSTERS.ConvexHall(i).numIons(1);
                y_ch = CLUSTERS.ConvexHall(i).numIons(2);
                if (cx==x_ch) && (cy==y_ch)
                    done = 1;
                    Eref = CLUSTERS.ConvexHall(i).Enthalpy;
                end
            end
            
            if ~done
                for i=1:length(CLUSTERS.ConvexHall) %numcomp
                    xch = CLUSTERS.ConvexHall(i).numIons(1);
                    ych = CLUSTERS.ConvexHall(i).numIons(2);
                    r(i) = sqrt((cx-xch)^2 + (cy-ych)^2);
                end
                [r_upor, n_upor] = sort(r);
                n_current = 2;
                while (~done)
                    n_current = n_current + 1;
                    for i = 1:n_current-2
                        for j = i+1:n_current-1
                            x1 = CLUSTERS.ConvexHall(n_upor(i)).numIons(1);
                            y1 = CLUSTERS.ConvexHall(n_upor(i)).numIons(2);
                            x2 = CLUSTERS.ConvexHall(n_upor(j)).numIons(1);
                            y2 = CLUSTERS.ConvexHall(n_upor(j)).numIons(2);
                            x3 = CLUSTERS.ConvexHall(n_upor(n_current)).numIons(1);
                            y3 = CLUSTERS.ConvexHall(n_upor(n_current)).numIons(2);
                            S = AreaTriangle(x1,y1,x2,y2,x3,y3);
                            S1 = AreaTriangle(cx,cy,x1,y1,x2,y2);
                            S2 = AreaTriangle(cx,cy,x1,y1,x3,y3);
                            S3 = AreaTriangle(cx,cy,x2,y2,x3,y3);
                            if (S~=0) && (S1+S2+S3 <= S+0.00000001)
                                E1 = CLUSTERS.ConvexHall(n_upor(i)).Enthalpy;
                                E2 = CLUSTERS.ConvexHall(n_upor(j)).Enthalpy;
                                E3 = CLUSTERS.ConvexHall(n_upor(n_current)).Enthalpy;
                                Eref=(det([x1 y1 E1;x2 y2 E2;x3 y3 E3])-cx*det([y1 E1 1;y2 E2 1;y3 E3 1])+cy*det([x1 E1 1;x2 E2 1;x3 E3 1]))/det([x1 y1 1;x2 y2 1;x3 y3 1]);
                                done = 1;
                            end
                            if done break; end
                        end
                        if done break; end
                    end
                end
                
            end
            n=n+1;
            xx(n) = cx + (1/ny) * (cy-numions(1,2));
            ff(n) = CLUSTERS.composition(cx-numions(1,1)+2,cy-numions(1,2)+2).bestEnthalpy - Eref;
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
                n = n + 1;
                x(n) = ix + (1/ny) * (iy-numions(1,2));
                F(n) = fitness(i);
            end
        end
        scatter(x,F,s,ccc(k,:),'filled','MarkerEdgeColor','k');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% variation operators %%%%%%%%%%%%%%%%%%%%%%%%%%
    min_ax = numions(1,1);
    max_ax = numions(2,1)+1;
    min_ay = 0;
    max_ay = 10; %CLUSTERS.fitness_range;
    xxx = max_ax - 0.5/(numions(2,2)-numions(1,2)+1); % - (max_ax-min_ax)/10;
    
    for i = 1:size(operators_text,1)
        yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
        scatter(xxx,yyy,60,ccc(i,:),'filled','MarkerEdgeColor','k');
        text(xxx+(max_ax-min_ax)/50, yyy, ['(' num2str(num_opera(i)) ') ' operators_text(i,:)]);
    end
    axis([min_ax,max_ax,min_ay,max_ay]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%% ticks & labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    n=0;
    max_ticks = 20;
    N1_ticks = fix(max_ticks/nx);
    delta_tick = fix(ny/N1_ticks) +1;
    if (delta_tick==0) delta_tick = 1; end
    for ix = numions(1,1):numions(2,1)
        for iy = numions(1,2):numions(2,2)
            n=n+1;
            x_tick(n) = ix + (1/ny) * (iy-numions(1,2));
            s_tick(n,:) = '     ';
            if (fix((iy-numions(1,2))/delta_tick)*delta_tick == iy-numions(1,2)) && (numions(2,2)+1-iy >= delta_tick)
                str_tick = [num2str(ix) ',' num2str(iy)];
                if length(num2str(ix))==1 str_tick = [' ' str_tick]; end
                if length(num2str(iy))==1 str_tick = [str_tick ' ']; end
                s_tick(n,:) = str_tick;
            end
        end
    end
    
    set(gca,'XTick', x_tick,'XTickLabel',s_tick);
    print(h,'-dtiff','-r120',[graphFolder '/Fitness_gen' num2str(POP_STRUC.generation) '.tif']);
    close(h);
catch
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Energy/numions vs numions %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% try
%     h = figure;
%     hold;
%     set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing
%     xlabel('numIons');
%     ylabel('Energy / numIons');
%     
%     % Plotting line connecting best clusters of all compositions in current generations %
%     Ebest = [];
%     Xbest = [];
%     n=0;
%     for ix = numions(1,1):numions(2,1)
%         for iy = numions(1,2):numions(2,2)
%             n=n+1;
%             Ebest(n) = CLUSTERS.composition(ix-numions(1,1)+2,iy-numions(1,2)+2).bestEnthalpy / (ix+iy);
%             Xbest(n) = ix + (1/ny) * (iy-numions(1,2));
%         end
%     end
%     plot(Xbest,Ebest,'LineWidth',1);
%     
%     %%%%%%%%%%%%%%%%% Clusters in current generation %%%%%%%%%%%%%%%%%%%%%
%     for k = 1:size(operators,1)
%         x=[]; EdivN=[]; n=0;
%         for i=1:length(POP_STRUC.POPULATION)
%             if ~isempty(strfind(POP_STRUC.POPULATION(i).howCome, operators(k,:)))
%                 ix = POP_STRUC.POPULATION(i).numIons(1);
%                 iy = POP_STRUC.POPULATION(i).numIons(2);
%                 n = n + 1;
%                 x(n) = ix + (1/ny) * (iy-numions(1,2));
%                 EdivN(n) = POP_STRUC.POPULATION(i).Enthalpies(end) / (ix+iy);
%             end
%         end
%         scatter(x,EdivN,s,ccc(k,:),'filled','MarkerEdgeColor','k');
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%% Plotting reference data %%%%%%%%%%%%%%%%%%%%%%
%     if (fid_ref~=-1)
%         Eref = [];
%         Xref = [];
%         n=0;
%         for ix = numions(1,1):numions(2,1)
%             for iy = numions(1,2):numions(2,2)
%                 if (E_LJ(ix,iy)~=0)
%                     n=n+1;
%                     Eref(n) = E_LJ(ix,iy) / (ix+iy);
%                     Xref(n) = ix + (1/ny) * (iy-numions(1,2));
%                     plot([Xref(n)-1/(2*ny), Xref(n)+1/(2*ny)], [Eref(n),Eref(n)],'k','LineWidth',1);
%                 end
%             end
%         end
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%% variation operators %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     min_ax = numions(1,1);
%     max_ax = numions(2,1)+1;
%     min_ay = min(Ebest)-(max(Ebest)-min(Ebest))/30;
%     max_ay = max(Ebest) + (max(Ebest)-min(Ebest))/2;
%     if (fid_ref~=-1) && (min(Eref)<min(Ebest))
%         min_ay = min(Eref)-(max(Ebest)-min(Eref))/30;
%         max_ay = max(Ebest) + (max(Ebest)-min(Eref))/2;
%     end
%     xxx = max_ax - 0.5/(numions(2,2)-numions(1,2)+1); % - (max_ax-min_ax)/10;
%     for i = 1:size(operators_text,1);
%         yyy = max_ay - (i-1)*(max_ay-min_ay)/30;
%         scatter(xxx,yyy,60,ccc(i,:),'filled','MarkerEdgeColor','k');
%         text(xxx+(max_ax-min_ax)/50, yyy, ['(' num2str(num_opera(i)) ') ' operators_text(i,:)]);
%     end
%     if (fid_ref~=-1)
%         plot([xxx-(max_ax-min_ax)/60, xxx+(max_ax-min_ax)/60], [yyy-(max_ay-min_ay)/30, yyy-(max_ay-min_ay)/30],'k','LineWidth',1);
%         text(xxx+(max_ax-min_ax)/40, yyy-(max_ay-min_ay)/30, 'Ref. data');
%     end
%     axis([min_ax,max_ax,min_ay,max_ay]);
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%% ticks & labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     n=0;
%     max_ticks = 20;
%     N1_ticks = fix(max_ticks/nx);
%     delta_tick = fix(ny/N1_ticks) +1;
%     if (delta_tick==0) delta_tick = 1; end
%     for ix = numions(1,1):numions(2,1)
%         for iy = numions(1,2):numions(2,2)
%             n=n+1;
%             x_tick(n) = ix + (1/ny) * (iy-numions(1,2));
%             s_tick(n,:) = '     ';
%             if (fix((iy-numions(1,2))/delta_tick)*delta_tick == iy-numions(1,2)) && (numions(2,2)+1-iy >= delta_tick)
%                 str_tick = [num2str(ix) ',' num2str(iy)];
%                 if length(num2str(ix))==1 str_tick = [' ' str_tick]; end
%                 if length(num2str(iy))==1 str_tick = [str_tick ' ']; end
%                 s_tick(n,:) = str_tick;
%             end
%         end
%     end
%     
%     set(gca,'XTick', x_tick,'XTickLabel',s_tick);
%     print(h,'-dtiff','-r120',[graphFolder '/Energy_div_numIons_gen' num2str(POP_STRUC.generation) '.tif']);
%     close(h);
% catch
% end

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
        n=0;
        for ix = numions(1,1):numions(2,1)
            for iy = numions(1,2):numions(2,2)
                n = n + 1;
                Erel(n) = CLUSTERS.composition(ix-numions(1,1)+2,iy-numions(1,2)+2).bestEnthalpy - E_LJ(ix,iy);
                Xrel(n) = ix + (1/ny) * (iy-numions(1,2));
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
                n = n + 1;
                Xmag(n) = ix + (1/ny) * (iy-numions(1,2));
                Emag(n) = POP_STRUC.POPULATION(i).Enthalpies(end) - E_LJ(ix,iy);
            end
        end
        scatter(Xmag,Emag,s,'k','filled','MarkerEdgeColor','k');
        
        %%%%%%%%%%%%%%%%%%%%% diapason of axis %%%%%%%%%%%%%%%%%%%%%%%%%%%
        min_ax = numions(1,1);
        max_ax = numions(2,1)+1;
        min_ay = -0.5;
        max_ay = 3;
        axis([min_ax,max_ax,min_ay,max_ay]);
        plot([min_ax,max_ax],[0,0],'k');
        
        %%%%%%%%%%%%%%%%%%%%%%% ticks & labels %%%%%%%%%%%%%%%%%%%%%%%%%%%
        n=0;
        max_ticks = 20;
        N1_ticks = fix(max_ticks/nx);
        delta_tick = fix(ny/N1_ticks) +1;
        if (delta_tick==0) delta_tick = 1; end
        for ix = numions(1,1):numions(2,1)
            for iy = numions(1,2):numions(2,2)
                n=n+1;
                x_tick(n) = ix + (1/ny) * (iy-numions(1,2));
                s_tick(n,:) = '     ';
                if (fix((iy-numions(1,2))/delta_tick)*delta_tick == iy-numions(1,2)) && (numions(2,2)+1-iy >= delta_tick)
                    str_tick = [num2str(ix) ',' num2str(iy)];
                    if length(num2str(ix))==1 str_tick = [' ' str_tick]; end
                    if length(num2str(iy))==1 str_tick = [str_tick ' ']; end
                    s_tick(n,:) = str_tick;
                end
            end
        end
        set(gca,'XTick', x_tick,'XTickLabel',s_tick);
        print(h,'-dtiff','-r120',[graphFolder '/Energy-refEnergy_gen' num2str(POP_STRUC.generation) '.tif']);
        close(h);
    catch
    end
end
