function WriteClusters(fitness)
global ORG_STRUC
global POP_STRUC
global CLUSTERS

resFolder = POP_STRUC.resFolder;
atomType= ORG_STRUC.atomType;
numions = ORG_STRUC.numIons;

graphFolder = [resFolder '/__graphsFitness_etal'];
if ~exist(graphFolder)
    mkdir(graphFolder)
end

bestFolder = [resFolder '/__thebestClusters'];
if ~exist(bestFolder)
    mkdir(bestFolder)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%------------------------------------ case of one-component clusters -------------------------------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(numions,2) == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% printing files best_clusters.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fpath  = [bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation)];
    fid  = fopen(fpath,  'wt');
    fprintf(fid,'composition  energy - all compositions\n');
    for i = 1:CLUSTERS.number_compositions
        fprintf(fid,'%9.0f  %9.3f\n',CLUSTERS.composition(i).numIons, CLUSTERS.composition(i).bestEnthalpy);
    end
    
    fprintf(fid,'\n');
    fprintf(fid,'composition  energy - convex hull\n');
    for i = 1:CLUSTERS.ConvexHall_numComp
        fprintf(fid,'%9.0f  %9.3f\n',CLUSTERS.ConvexHall(i).numIons, CLUSTERS.ConvexHall(i).Enthalpy);
    end
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% printing files best_clusters_POSCARS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = numions(1,1) : numions(2,1)
        ii = i-numions(1,1)+1;
        coor     = CLUSTERS.composition(ii).COORDINATES;
        numIons1 = CLUSTERS.composition(ii).numIons;
        lattice  = CLUSTERS.composition(ii).LATTICE;
        symg     = CLUSTERS.composition(ii).symg;
        order    = CLUSTERS.composition(ii).order;
        count    = CLUSTERS.composition(ii).Number;
        Write_POSCAR(atomType, count, symg, numIons1, lattice, coor);
        unixCmd([' cat POSCAR       >>' bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation) '_POSCARS']);
    end
    
    plot_graphs_1comp(graphFolder, fitness);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------case of two-component (binary) clusters-----------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(numions,2) == 2
    
    %%%%%%%%%%%%% printing files best_clusters.txt %%%%%%%%%%%%%%%%%%%%%%%
    fpath  = [bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation)];
    fid  = fopen(fpath,  'wt');
    
    %printing energies of the best clusters
    fprintf(fid,'        |');
    for i = numions(1,2) : numions(2,2)
        fprintf(fid,'    %2.0f     ',i);
    end
    fprintf(fid,'\n');
    for i = numions(1,2) : numions(2,2) + 1
        fprintf(fid,'----------');
    end
    fprintf(fid,'\n');
    
    for i = numions(1,1) : numions(2,1)
        fprintf(fid,'  %2.0f    |',i);
        for j = numions(1,2) : numions(2,2)
            ii = i-numions(1,1)+2;
            jj = j-numions(1,2)+2;
            fprintf(fid,' %9.4f ',CLUSTERS.composition(ii,jj).bestEnthalpy);
        end
        fprintf(fid,'\n');
    end
    
    %printing magic clusters
    fprintf(fid,'\n');
    fprintf(fid,'    |');
    for i = numions(1,2) : numions(2,2)
        fprintf(fid,' %2.0f ',i);
    end
    fprintf(fid,'\n');
    for i = numions(1,2) : numions(2,2) + 1
        fprintf(fid,'----');
    end
    fprintf(fid,'\n');
    
    for i = numions(1,1) : numions(2,1)
        fprintf(fid,' %2.0f |',i);
        for j = numions(1,2) : numions(2,2)
            ii = i-numions(1,1)+2;
            jj = j-numions(1,2)+2;
            d2Ex = CLUSTERS.composition(ii+1,jj).bestEnthalpy+CLUSTERS.composition(ii-1,jj).bestEnthalpy-2*CLUSTERS.composition(ii,jj).bestEnthalpy;
            d2Ey = CLUSTERS.composition(ii,jj+1).bestEnthalpy+CLUSTERS.composition(ii,jj-1).bestEnthalpy-2*CLUSTERS.composition(ii,jj).bestEnthalpy;
            if (d2Ex >= 0) && (d2Ey >= 0)
                fprintf(fid,'  m ');
            else
                fprintf(fid,'    ');
            end
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    %%%%%%%%%%%% printing files best_clusters_POSCARS.txt %%%%%%%%%%%%%%%%
    
    for i = numions(1,1) : numions(2,1)
        for j = numions(1,2) : numions(2,2)
            ii = i-numions(1,1)+2;
            jj = j-numions(1,2)+2;
            coor     = CLUSTERS.composition(ii,jj).COORDINATES;
            numIons1 = CLUSTERS.composition(ii,jj).numIons;
            lattice  = CLUSTERS.composition(ii,jj).LATTICE;
            symg     = CLUSTERS.composition(ii,jj).symg;
            order    = CLUSTERS.composition(ii,jj).order;
            count    = CLUSTERS.composition(ii,jj).Number;
            Write_POSCAR(atomType, count, symg, numIons1, lattice, coor);
            unixCmd([' cat POSCAR       >>' bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation) '_POSCARS']);
        end
    end
    
    plot_graphs_2comp(graphFolder, fitness);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%--------------case of three-component (binary) clusters---------------%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(numions,2) == 3
    el1 = megaDoof(ORG_STRUC.atomType(1));
    el2 = megaDoof(ORG_STRUC.atomType(2));
    el3 = megaDoof(ORG_STRUC.atomType(3));

    %%%%%%%%%%%%% printing files best_clusters.txt %%%%%%%%%%%%%%%%%%%%%%%
    fpath  = [bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation)];
    fid  = fopen(fpath,  'wt');

    %printing energies of the best clusters
    fprintf(fid,'Best energies:\n');
    for i = numions(1,1) : numions(2,1)
        ii = i-numions(1,1)+2;
        
        fprintf(fid,['  ' el1 num2str(i) '   |']);
        for k = numions(1,3) : numions(2,3)
            fprintf(fid,['    ' el3 num2str(k) '    '],i);
        end
        fprintf(fid,'\n');
        for k = numions(1,3) : numions(2,3) + 1
            fprintf(fid,'----------');
        end
        fprintf(fid,'\n');
        
        for j = numions(1,2) : numions(2,2)
            fprintf(fid,['  ' el2 num2str(j) '  |']);
            for k = numions(1,3) : numions(2,3)
                jj = j-numions(1,2)+2;
                kk = k-numions(1,3)+2;
                fprintf(fid,' %9.4f ',CLUSTERS.composition(ii,jj,kk).bestEnthalpy);
            end
            fprintf(fid,'\n');
        end
        
        fprintf(fid,'\n');
    end
    
    %printing magic clusters
    fprintf(fid,'Magic clusters:\n');
    for i = numions(1,1) : numions(2,1)
        ii = i-numions(1,1)+2;
        
        fprintf(fid,[' ' el1 num2str(i) ' |']);
        for k = numions(1,3) : numions(2,3)
            fprintf(fid,[' ' el3 num2str(k) '  '],i);
        end
        fprintf(fid,'\n');
        for k = numions(1,3) : numions(2,3) + 1
            fprintf(fid,'-----');
        end
        fprintf(fid,'\n');
        
        for j = numions(1,2) : numions(2,2)
            fprintf(fid,[' ' el2 num2str(j) ' |']);
            for k = numions(1,3) : numions(2,3)
                jj = j-numions(1,2)+2;
                kk = k-numions(1,3)+2;
                
                d2Ex = CLUSTERS.composition(ii+1,jj,kk).bestEnthalpy+CLUSTERS.composition(ii-1,jj,kk).bestEnthalpy-2*CLUSTERS.composition(ii,jj,kk).bestEnthalpy;
                d2Ey = CLUSTERS.composition(ii,jj+1,kk).bestEnthalpy+CLUSTERS.composition(ii,jj-1,kk).bestEnthalpy-2*CLUSTERS.composition(ii,jj,kk).bestEnthalpy;
                d2Ez = CLUSTERS.composition(ii,jj,kk+1).bestEnthalpy+CLUSTERS.composition(ii,jj,kk-1).bestEnthalpy-2*CLUSTERS.composition(ii,jj,kk).bestEnthalpy;
                if (d2Ex >= 0) && (d2Ey >= 0) && (d2Ez >= 0)
                    fprintf(fid,'  m  ');
                else
                    fprintf(fid,'     ');
                end
            end
            fprintf(fid,'\n');
        end

        fprintf(fid,'\n');
    end
    
    fclose(fid);
    
    %%%%%%%%%%%% printing files best_clusters_POSCARS.txt %%%%%%%%%%%%%%%%
    
    for i = numions(1,1) : numions(2,1)
    for j = numions(1,2) : numions(2,2)
    for k = numions(1,3) : numions(2,3)
        ii = i-numions(1,1)+2;
        jj = j-numions(1,2)+2;
        kk = k-numions(1,3)+2;
        coor     = CLUSTERS.composition(ii,jj,kk).COORDINATES;
        numIons1 = CLUSTERS.composition(ii,jj,kk).numIons;
        lattice  = CLUSTERS.composition(ii,jj,kk).LATTICE;
        symg     = CLUSTERS.composition(ii,jj,kk).symg;
        order    = CLUSTERS.composition(ii,jj,kk).order;
        count    = CLUSTERS.composition(ii,jj,kk).Number;
        Write_POSCAR(atomType, count, symg, numIons1, lattice, coor);
        unixCmd([' cat POSCAR       >>' bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation) '_POSCARS']);
    end
    end
    end
    
    plot_graphs_3comp(graphFolder, fitness);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%creating files best_clusters and best_clusters_POSCARS in "results" folder
unixCmd(['cat /dev/null > ' resFolder '/!best_clusters']);
unixCmd(['cat /dev/null > ' resFolder '/!best_clusters_POSCARS']);
unixCmd([' cat ' bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation) '       >>' resFolder '/!best_clusters']);
unixCmd([' cat ' bestFolder '/best_clusters_gen' num2str(POP_STRUC.generation) '_POSCARS       >>' resFolder '/!best_clusters_POSCARS']);
