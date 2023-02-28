function plot_all_isomers()
global ORG_STRUC
global POP_STRUC
global USPEX_STRUC
global CLUSTERS

graphFolder = [POP_STRUC.resFolder '/___graphs_Energy_Fitness___'];
numions = ORG_STRUC.numIons;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------ case of one-component clusters -------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(numions,2) == 2

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
            E_RefData(i,j) = fscanf(fid_ref,'%f',1);
        end
        ss = fgets(fid_ref);
    end
    fclose(fid_ref);

    try
        %%%%%%%%%%%%%%%%%%%%%%%% plot the graph %%%%%%%%%%%%%%%%%%%%%%%%%%
        h = figure;
        hold;
        set(gcf,'Visible','off'); % Use this switcher to prevent Matlab foregroung printing
        xlabel('numIons');
        ylabel('Energy - reference Energy');
        
        %%%%%%%%%%%%%%%% diapason of axis; line y=0 %%%%%%%%%%%%%%%%%%%%%%
        min_ax = numions(1,1) - 1/(2*ny);
        max_ax = numions(2,1) + 1 + 1/(2*ny);
        min_ay = -1;
        max_ay = 5;
        axis([min_ax,max_ax,min_ay,max_ay]);
        plot([min_ax,max_ax],[0,0],'k:','LineWidth',0.5);
        
        for ix = numions(1,1):numions(2,1)
            for iy = numions(1,2):numions(2,2)
                Xrel = ix + (1/ny) * (iy-numions(1,2));
                plot([Xrel + 1/(2*ny), Xrel + 1/(2*ny)], [min_ay,max_ay],'k:','LineWidth',0.5);
            end
        end
                
        % Plotting all isomers from CLUSTER.Structure database %
        for i = 1 : length(CLUSTERS.Structure)
            E = CLUSTERS.Structure(i).enthalpy;
            curIons = CLUSTERS.Structure(i).numIons;
            ix = curIons(1);
            iy = curIons(2);
            Erel = E - E_RefData(ix,iy); %energy - ref energy
            Xrel = ix + (1/ny) * (iy-numions(1,2));
            plot([Xrel - 1/(2*ny), Xrel], [Erel,Erel],'k','LineWidth',1);
        end
        
        % Plotting all isomers from USPEX_STRUC
        N = length(USPEX_STRUC.POPULATION);
        
        for ix = numions(1,1) : numions(2,1)
            for iy = numions(1,2) : numions(2,2)
                curions = [ix, iy];
                enth = [];
                number = []; % numner of current structure in USPEX_STRUC (array)
                nn = 0;      % the total number of structures in USPEX_STRUC with current composition
                for i = 1 : N
                    if (USPEX_STRUC.POPULATION(i).numIons == curions)
                        nn = nn+1;
                        number(nn) = i;
                        enth(nn) = USPEX_STRUC.POPULATION(i).Enthalpies(end);
                    end
                end
                [nothing, ranking] = sort(enth);
                
                % Remove all duplicates:
                %--------------------------------------------------------------------------
                GoodList = [];
                item = 1;
                GoodList(1) = 1; % Good Structure
                for i = 2 : nn
                    ii = number(ranking(i)); %number in USPEX_STRUC
                    if USPEX_STRUC.POPULATION(ii).ToCount > 0
                        same = 0;
                        for j= 1 : i-1
                            jj = number(ranking(j)); %number in USPEX_STRUC
                            if USPEX_STRUC.POPULATION(jj).ToCount > 0
                                if SameStructure(ii,jj,USPEX_STRUC)
                                    USPEX_STRUC.POPULATION(ii).ToCount = 0;
                                    same = 1;
                                    break;
                                end
                            end
                        end
                        if same == 0
                            item = item + 1;
                            GoodList(item) = i;
                        end
                    end
                end
                
                Xrel = ix + (1/ny) * (iy-numions(1,2));
                for ii = 1:item
                    Erel = enth(ranking(GoodList(ii))) - E_RefData(ix,iy); %energy - ref energy
                    plot([Xrel, Xrel + 1/(2*ny)], [Erel,Erel],'k','LineWidth',1);
                end
            end
        end
       
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
        print(h,'-dtiff','-r120',[graphFolder '/All_isomers_gen' num2str(POP_STRUC.generation) '.tif']);
    catch
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------case of two-component (binary) clusters-----------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(numions,2) == 1
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
    
        try
            %%%%%%%%%%%%%%%%%%%%%%%% plot the graph %%%%%%%%%%%%%%%%%%%%%%
            h = figure;
            hold;
            set(gcf,'Visible','off'); % Use this switcher to prevent Matlab foregroung printing
            xlabel('numIons');
            ylabel('Energy - reference Energy');
            
            %%%%%%%%%%%%%%%%%%%%%%% diapason of axis %%%%%%%%%%%%%%%%%%%%%
            min_ax = numions(1,1) - 0.5;
            max_ax = numions(2,1) + 0.5;
            min_ay = -1;
            max_ay = 5;
            axis([min_ax,max_ax,min_ay,max_ay]);
            plot([min_ax,max_ax],[0,0],'k:','LineWidth',0.5);
            
            for ix = numions(1,1) : numions(2,1)
                plot([ix+0.5, ix+0.5], [min_ay,max_ay],'k:','LineWidth',0.5);
            end
            
            % Plotting all isomers from CLUSTER.Structure database %
            for i = 1 : length(CLUSTERS.Structure)
                ix = CLUSTERS.Structure(i).numIons;
                Erel = CLUSTERS.Structure(i).enthalpy - E_RefData(ix); %energy - ref energy
                plot([ix - 0.5, ix], [Erel,Erel],'k','LineWidth',1);
            end
            
            % Plotting all isomers from USPEX_STRUC
            N = length(USPEX_STRUC.POPULATION);
            
            for ix = numions(1,1) : numions(2,1)
                enth = [];
                number = []; % numner of current structure in USPEX_STRUC (array)
                nn = 0;      % the total number of structures in USPEX_STRUC with current composition
                for i = 1 : N
                    if (USPEX_STRUC.POPULATION(i).numIons == ix)
                        nn = nn+1;
                        number(nn) = i;
                        enth(nn) = USPEX_STRUC.POPULATION(i).Enthalpies(end);
                    end
                end
                [nothing, ranking] = sort(enth);
                
                % Remove all duplicates:
                %--------------------------------------------------------------------------
                GoodList = [];
                item = 1;
                GoodList(1) = 1; % Good Structure
                for i = 2 : nn
                    ii = number(ranking(i)); %number in USPEX_STRUC
                    if USPEX_STRUC.POPULATION(ii).ToCount > 0
                        same = 0;
                        for j= 1 : i-1
                            jj = number(ranking(j)); %number in USPEX_STRUC
                            if USPEX_STRUC.POPULATION(jj).ToCount > 0
                                if SameStructure(ii,jj,USPEX_STRUC)
                                    USPEX_STRUC.POPULATION(ii).ToCount = 0;
                                    same = 1;
                                    break;
                                end
                            end
                        end
                        if same == 0
                            item = item + 1;
                            GoodList(item) = i;
                        end
                    end
                end
                
                for ii = 1:item
                    Erel = enth(ranking(GoodList(ii))) - E_RefData(ix); %energy - ref energy
                    plot([ix, ix + 0.5], [Erel,Erel],'k','LineWidth',1);
                end
                print(h,'-dtiff','-r120',[graphFolder '/All_isomers_gen' num2str(POP_STRUC.generation) '.tif']);
            end
        catch
        end
    end    
end