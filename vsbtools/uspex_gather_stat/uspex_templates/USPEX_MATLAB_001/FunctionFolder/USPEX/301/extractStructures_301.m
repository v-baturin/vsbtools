function extractStructures_301(max_num, weight)
global POOL_STRUC
global USPEX_STRUC
global ORG_STRUC
% This function creates a list of structures that are closer than fit_max to the best
% Identical structures are kicked out by enthalpy and order parameter
% Last updated by Qiang Zhu (2013/12/27)
% read enthalpies and compositions
Ntype = length(USPEX_STRUC.POPULATION(1).numIons);
N     = length(USPEX_STRUC.POPULATION);
enth    = zeros(N,1);
fitness = zeros(N,1);
comp    = zeros(N,Ntype);

for i=1:N
    comp(i,:)  = USPEX_STRUC.POPULATION(i).numIons;
    fitness(i) = USPEX_STRUC.POPULATION(i).Fitness;
    enth(i)    = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(comp(i,:));
end
if ORG_STRUC.paretoRanking == 0
    [nothing, ranking] = sort(fitness);
else
%%%% Here First, we Just rank systems with respect to convex hull
%%%% Then discard bad ones and keep BestFrac of them, And then 
%%%% Rank remained structures using Pareto Method for more than one Property
%%%% This is suggested by Artem, Implemented by Zahed ZX.
    convex = distance_convex_hull(ORG_STRUC.numIons, ORG_STRUC.homePath);
    [nothing, RANK] = sort(convex);
    N_pro = round(length(convex) * ORG_STRUC.bestFrac);
    badRank = RANK(N_pro+1:end);

    [ranking, ParetoFront, Properties, MSG] = ParetoRankingAll();   
    for rr = 1 : length(badRank)
        for a = length(ranking) : -1 : 1
            if ranking(a) == badRank(rr)
		ranking(a) = [];
	    end
	end
    end
end
%--------------------------------------------------------------------------
% First step: fine filtering:
N = length(ranking);
for i = 2 : N
    if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
        for j= 1 : i-1
            if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
                if SameStructure(ranking(i), ranking(j), USPEX_STRUC)
                    USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
                    break;
                end
            end
        end
    end
end
%--------------------------------------------------------------------------
% Second step: rough filtering using A_order:
GoodList = zeros(max_num, 1); % Good Structure
item = 1;
GoodList(1) = 1; % Good Structure
for i = 2 : N
    if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
        same = 0;
        for j= 1 : i-1
            if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
                if SameStructure_order(ranking(i), ranking(j), USPEX_STRUC)
                    USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
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
    
    if item == max_num
        break;
    end
end
%--------------------------------------------------------------------------

Good_num = item; %less than max_num

% output all accepted structures
fp1 = fopen('goodStructures', 'w');
fp2 = fopen('goodStructures_POSCARS', 'w');
fprintf(fp1,'  ID   Compositions    Enthalpies    Volumes    fitness  SYMM\n');
fprintf(fp1,'                       (eV/atom)    (A^3/atom)    ()         \n');
for count = 1 : Good_num
    j = ranking(GoodList(count));
    symg      = USPEX_STRUC.POPULATION(j).symg;
    num       = USPEX_STRUC.POPULATION(j).numIons;
    lat       = USPEX_STRUC.POPULATION(j).LATTICE;
    coords    = USPEX_STRUC.POPULATION(j).COORDINATES;
    order     = USPEX_STRUC.POPULATION(j).order;
    numBlocks = USPEX_STRUC.POPULATION(j).numBlocks;
    fing      = USPEX_STRUC.POPULATION(j).FINGERPRINT;
    atomType_str = USPEX_STRUC.SYSTEM.atomType_symbol;
    volume = det(lat)/sum(num);
    
    composition = sprintf('%3d',num);
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
        composition=[composition,blanks(shift(length(num)))];
    end
    
    fprintf(fp1,'%4d  [%11s]   %9.4f   %9.4f   %9.4f  %4d\n', j, composition, enth(j), volume, fitness(j), symg);
    % POSCARS
    fprintf(fp2,'Structure_%-4d Sym.group: %4d\n', j, symg);
    fprintf(fp2, '1.0000\n');
    
    for latticeLoop = 1 : 3
        fprintf(fp2, '%12.6f %12.6f %12.6f\n', lat(latticeLoop,:));
    end
    for i = 1 : length(atomType_str)
        fprintf(fp2, '%4s ', atomType_str{i});
    end
    fprintf(fp2, '\n');
    for i=1:length(num)
        fprintf(fp2, '%4d ', num(i));
    end
    fprintf(fp2, '\n');
    fprintf(fp2, 'Direct\n');
    for coordLoop = 1 : sum(num)
        fprintf(fp2, '%12.6f %12.6f %12.6f\n', coords(coordLoop,:));
    end
    
    POOL_STRUC.POPULATION(count).LATTICE = lat;
    POOL_STRUC.POPULATION(count).COORDINATES = coords;
    POOL_STRUC.POPULATION(count).numIons = num;
    POOL_STRUC.POPULATION(count).numBlocks = numBlocks;
    POOL_STRUC.POPULATION(count).order = order;
    POOL_STRUC.POPULATION(count).enthalpy = enth(j);
    POOL_STRUC.POPULATION(count).Number = j;
    POOL_STRUC.POPULATION(count).FINGERPRINT = fing;
end
if ORG_STRUC.paretoRanking ~= 0					    
    for a = 1 : Good_num                                            
	ParetoF(a,1) = ParetoFront(GoodList(a));                    
    end                                                             
    %%% This part is to put structures from upper Pareto front      
    %%% to lower one when one Pareto front doesn't have Structures  
    CheckPoint = 1;                                                 
    for a = 1 : length(ParetoF)                                 
        if ParetoF(a) == CheckPoint                             
        elseif ParetoF(a) == CheckPoint + 1                     
            CheckPoint = CheckPoint + 1;                        
        else                                                    
    	diff = ParetoF(a) - (CheckPoint + 1);               
    	ParetoF(a:end) = ParetoF(a:end) - diff;             
    	CheckPoint = CheckPoint + 1;
        end                                                     
    end                                                         
    %%------------------------------------------------------------
    ParetoFrontier = zeros(ParetoF(end),1);
    for a = 1 : length(ParetoFrontier)   %% equal to number of Pareto front lines
        for b = 1 : length(ParetoF)  %% equal to number of structures
            if ParetoF(b) == a
                ParetoFrontier(a) = ParetoFrontier(a) + 1;
            end
        end
    end
    POOL_STRUC.paretoFront = ParetoFrontier;
end
fclose(fp1);
fclose(fp2);
