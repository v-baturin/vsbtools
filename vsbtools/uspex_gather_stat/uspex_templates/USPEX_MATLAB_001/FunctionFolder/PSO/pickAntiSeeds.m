function pickAntiSeeds()

% USPEX Version 7.4.3
% Change: added

global ANTISEEDS
global ORG_STRUC

ANTISEEDS = struct('FINGERPRINT',{},'Sigma',{},'Max',{});
ANTISEEDS(1).FINGERPRINT = [];
ANTISEEDS(1).Sigma = [];
ANTISEEDS(1).Max = [];

cd AntiSeeds

try
    
    %%%%%%%%%%% get atom positions %%%%%%%%%%%%%%%%%%
    
    % reading compositions
    if ORG_STRUC.varcomp == 1
        [fid,message] = fopen('compositions','r+');
        compos = fscanf(fid, '%d', [length(ORG_STRUC.atomType) inf]);
        compos = compos';
        status = fclose(fid);
    else
        compos = ORG_STRUC.numIons;
    end
    
    % reading POSCARS
    
    [fid,message] = fopen('POSCARS');
    loops_seed = 0;
    while 1
        try
            
            loops_seed = loops_seed + 1;
            
            tmp = fgetl(fid); % system description
            scale_factor = fgetl(fid); % 1 = numbers in angstrems
            optlattice = fscanf(fid,'%g',[3,3]); % lattice vectors
            optlattice = optlattice';
            
            tmp = fgetl(fid);  % we don't need this line
            ntyp = fgetl(fid); % types of atoms aka atomic numbers
            ntyp = str2num(ntyp);
            natom = sum(ntyp); % number of atoms
            tmp = fgetl(fid);
            candidate_pop = fscanf(fid,'%g',[3,natom]);
            candidate_pop = candidate_pop';
            tmp = fgetl(fid);
            
            tic
            if ORG_STRUC.varcomp == 1
                numIons = compos(loops_seed,:);
            else
                numIons = ORG_STRUC.numIons;
            end
            natom = sum(numIons);
            if ORG_STRUC.varcomp == 1
                [A, Ni, V, Nfull, W, dist_matrix, typ_i, typ_j, species] = makeMatrices(ORG_STRUC.RmaxFing, optlattice, natom, candidate_pop, 1, natom, 1);
            else
                [A, Ni, V, Nfull, W, dist_matrix, typ_i, typ_j, species] = makeMatrices(ORG_STRUC.RmaxFing, optlattice, numIons, candidate_pop, 1, natom, ORG_STRUC.atomType);
            end
            [order, ANTISEEDS(loops_seed).FINGERPRINT, atom_fing] = fingerprint_calc(ORG_STRUC.RmaxFing, Ni, natom, Nfull, A, W, V, dist_matrix, ORG_STRUC.sigmaFing, ORG_STRUC.deltaFing, species, typ_i, typ_j);
            fingCalcTime = toc    % shows calculation time for each fingerprint
            
        catch
            break
        end
    end
    status = fclose(fid);
    
catch
end

cd (ORG_STRUC.homePath)
