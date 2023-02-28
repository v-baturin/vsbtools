function coordinating_number_analisis()
%this function perform analisis of coordinating numbers of best structures in last processed generation
%it stores minimum and maximum coordinating number for each atom type and each composition
%Note: actual coordinating numbers evaluation performed in python script USPEX/src/calculate_coordinating_numbers.py
%developer - Bushlanov Pavel

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

howManyProliferate = min( [round((length(POP_STRUC.ranking))*ORG_STRUC.bestFrac),...
                                  length(POP_STRUC.ranking)-POP_STRUC.bad_rank ] );

atomType  = ORG_STRUC.atomType;


R_val = zeros(1,length(atomType));
for i = 1 : length(atomType)
    s = covalentRadius(ceil(atomType(i)));
    R_val(i) = str2num(s);
end
str_R_val = num2str(R_val);

compositions_str = cell(0);   %list of all compositions of best structures in last processed generation
compositions_cnbt = cell(0);  %list of all coordinating numbers of best structures in last processed generation
type_number = size(atomType,2);
coordinating_numbers_by_type_min = cell(0); 
coordinating_numbers_by_type_max = cell(0);

for i = 1:howManyProliferate %for all the good Structures
    lattice = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).LATTICE;
    coordinates = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).COORDINATES;
    numIons = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).numIons;
    if sum(numIons) ~= size(coordinates,1)
	CoeF = size(coordinates,1) / sum(numIons);
	numIons = CoeF * numIons;
    end
    str_numIons = num2str(numIons);
    ind_LOCAL = find(strcmp(compositions_str, str_numIons));
    if isempty(ind_LOCAL)
        ind_LOCAL = size(compositions_str,1)+1;
        compositions_str{ind_LOCAL,1} = str_numIons;
        for j = 1:type_number
            compositions_cnbt{ind_LOCAL,j} = zeros (0, numIons(j));
        end
    end
    str_lattice = sprintf('% 10.5f', reshape(lattice.',1,[]));
    str_coordinates = sprintf('% 10.5f', reshape(coordinates.',1,[]));
    str_type_number = num2str(type_number);

    coordinating_numbers = python_uspex([ORG_STRUC.USPEXPath '/FunctionFolder/USPEX/src/calculate_coordinating_numbers.py'],str_type_number,str_numIons,str_R_val,str_lattice,str_coordinates,1);
    if size(coordinating_numbers,2) == 1
        return
    end
    atom_type_ranges = cumsum(horzcat([1],numIons));
    for j = 1:type_number
        compositions_cnbt{ind_LOCAL,j} = vertcat(compositions_cnbt{ind_LOCAL,j},coordinating_numbers(atom_type_ranges(j):atom_type_ranges(j+1)-1));
    end
end

%USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS is global structure containing information on 
%all coordinating numbers obtained in simulations
for ind_LOCAL = 1:size(compositions_str,1)
    for i = 1:type_number
        A = reshape(compositions_cnbt{ind_LOCAL,i},[],1);
        if ~isempty(A)
            coordinating_numbers_by_type_min{i} = min(A);
            coordinating_numbers_by_type_max{i} = max(A);
        else
            coordinating_numbers_by_type_min{i} = 3;
            coordinating_numbers_by_type_max{i} = 12;
        end
        if (ceil(coordinating_numbers_by_type_min{i})==2)&&(ceil(coordinating_numbers_by_type_max{i})==2)
            coordinating_numbers_by_type_min{i} = 3;
            coordinating_numbers_by_type_max{i} = 12;
        end
    end
    cnbtmin = cell2mat(coordinating_numbers_by_type_min);
    cnbtmax = cell2mat(coordinating_numbers_by_type_max);
    cmps_str = cell(0);
    for j = 1:size(USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS,1)
        cmps_str{j} = num2str(USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS(j).COMPOSITION);
    end
    ind_USPEX = find(strcmp(cmps_str, compositions_str{ind_LOCAL}));
    if isempty(ind_USPEX)
        ind_USPEX = size(USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS,2)+1;
        USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS(ind_USPEX).COMPOSITION = str2num(compositions_str{ind_LOCAL,1});
    end
    USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS(ind_USPEX).CN = vertcat(cnbtmin,cnbtmax);
end

%USPEX_STRUC.EXT_RandTop.CN_grid is global structure containing
%interpolated grid of coordinationg numbers for all compositions
X = vertcat(USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS(:).COMPOSITION)';
if isempty(ORG_STRUC.maxAt)
    grid_sz = sum(ORG_STRUC.numIons)+1;
else
    grid_sz = ORG_STRUC.maxAt+1;
end
if size(X,2)>3
    for i = 1:type_number
        V1 = cell2mat(cellfun(@(x) x(1,i),{USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS.CN},'UniformOutput',false));
        V2 = cell2mat(cellfun(@(x) x(2,i),{USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS.CN},'UniformOutput',false));
        switch type_number
            case 1
                X1q = 0:1:grid_sz;
                USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = ceil(interp1(X(1,:),V1,X1q,'nearest'));
                USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = ceil(interp1(X(1,:),V2,X1q,'nearest'));
            case 2
                [X1q,X2q] = meshgrid(0:1:grid_sz,0:1:grid_sz);
                USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = ceil(griddata(X(1,:),X(2,:),V1,X1q,X2q,'nearest'));
                USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = ceil(griddata(X(1,:),X(2,:),V2,X1q,X2q,'nearest'));
            case 3
                [X1q,X2q,X3q] = meshgrid(0:1:grid_sz,0:1:grid_sz,0:1:grid_sz);
                USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = ceil(griddata(X(1,:),X(2,:),X(3,:),V1,X1q,X2q,X3q,'nearest'));
                USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = ceil(griddata(X(1,:),X(2,:),X(3,:),V2,X1q,X2q,X3q,'nearest'));
            case 4
                [X1q,X2q,X3q,X4q] = meshgrid(0:1:grid_sz,0:1:grid_sz,0:1:grid_sz,0:1:grid_sz);
                USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = ceil(griddata(X(1,:),X(2,:),X(3,:),X(4,:),V1,X1q,X2q,X3q,X4q,'nearest'));
                USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = ceil(griddata(X(1,:),X(2,:),X(3,:),X(4,:),V2,X1q,X2q,X3q,X4q,'nearest'));
            case 5
                [X1q,X2q,X3q,X4q,X5q] = meshgrid(0:1:grid_sz,0:1:grid_sz,0:1:grid_sz,0:1:grid_sz,0:1:grid_sz);
                USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = ceil(griddata(X(1,:),X(2,:),X(3,:),X(4,:),X(5,:),V1,X1q,X2q,X3q,X4q,X5q,'nearest'));
                USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = ceil(griddata(X(1,:),X(2,:),X(3,:),X(4,:),X(5,:),V2,X1q,X2q,X3q,X4q,X5q,'nearest'));
            otherwise
                a_size = repmat([grid_sz],1,type_number);
                a_grid = ceil(repmat([3],a_size));
                USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = a_grid;
                a_grid = ceil(repmat([12],a_size));
                USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = a_grid;
        end
        if isempty(USPEX_STRUC.EXT_RandTop.CN_grid{1,i})
            a_size = repmat([grid_sz],1,type_number);
            a_grid = ceil(repmat([3],a_size));
            USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = a_grid;
        end
        if isempty(USPEX_STRUC.EXT_RandTop.CN_grid{2,i})
            a_size = repmat([grid_sz],1,type_number);
            a_grid = ceil(repmat([12],a_size));
            USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = a_grid;
        end
    end
else
    for i = 1:type_number
        a_size = repmat([grid_sz],1,type_number);
        a_grid = ceil(repmat([USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS(1).CN(1,i)],a_size));
        USPEX_STRUC.EXT_RandTop.CN_grid{1,i} = a_grid;
        a_grid = ceil(repmat([USPEX_STRUC.EXT_RandTop.COORDINATING_NUMBERS(1).CN(2,i)],a_size));
        USPEX_STRUC.EXT_RandTop.CN_grid{2,i} = a_grid;
    end
end
