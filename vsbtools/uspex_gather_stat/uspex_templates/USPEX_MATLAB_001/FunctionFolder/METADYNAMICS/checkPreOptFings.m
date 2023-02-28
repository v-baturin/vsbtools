function checkPreOptFings()

global ORG_STRUC
global POP_STRUC

try
    % displays the distances between structures before optimization
    unixCmd(['echo ----- generation' num2str(POP_STRUC.generation) ' ----- >> preOptDistances.dat']);
    distances = zeros(length(POP_STRUC.POPULATION), length(POP_STRUC.POPULATION));
    %natom = sum(POP_STRUC.POPULATION(1).numIons);
    %numIons = POP_STRUC.POPULATION(1).numIons;
    
    tic
    
    for i = 1 : length(POP_STRUC.POPULATION)
        lat = POP_STRUC.POPULATION(i).LATTICE;
        coor = POP_STRUC.POPULATION(i).COORDINATES;
        numIons = POP_STRUC.POPULATION(i).numIons;
        [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(lat, coor, numIons, ORG_STRUC.atomType);
        [order0, fingerprint0, atom_fing0] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
        for j = i+1 : length(POP_STRUC.POPULATION)
            lat1 = POP_STRUC.POPULATION(j).LATTICE;
            coor1 = POP_STRUC.POPULATION(j).COORDINATES;
            natom = sum(POP_STRUC.POPULATION(j).numIons);
            numIons = POP_STRUC.POPULATION(j).numIons;
            [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(lat1, coor1, numIons, ORG_STRUC.atomType);
            [order1, fingerprint1, atom_fing1] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons);
            d(i,j) = cosineDistance(fingerprint0, fingerprint1, ORG_STRUC.weight);
            d(j,i) = d(i,j);
        end
    end
    
    t = toc;
    disp(['Distance matrix calculated in ' num2str(t) 'sec']);
    
    if ~isempty(d)
        unixCmd(['echo ----- generation' num2str(POP_STRUC.generation) ' ----- >> ' ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/non_optimized_POSCARs']);
        for i = 1 : size(d,1)
            unixCmd([' echo ' num2str(d(i,:)) ' >> preOptDistances.dat']);
            %  writeOUT_POSCAR(i, 0);
            
            lattice = POP_STRUC.POPULATION(Ind_No).LATTICE;
            coord   = POP_STRUC.POPULATION(Ind_No).COORDINATES;
            numIons = POP_STRUC.POPULATION(Ind_No).numIons;
            count = POP_STRUC.POPULATION(Ind_No).Number;
            
            Write_POSCAR(ORG_STRUC.atomType, count, symg, numIons, lattice, coord);
            
            unixCmd([' cat POSCAR  >> ' ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/non_optimized_POSCARs']);
        end
    end
    
catch
    disp('Error, during calculation of the distance matrix before optimization! Possible reason: number of coordinates is not equal to sum(numIons) for some structure.');
end
