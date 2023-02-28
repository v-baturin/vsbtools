function [all_freqs, all_eigvectors, supercells_sizes] = calcSoftModes_varcomp(lat, coords, numIons, maxAtoms, maxIncrease)

% calculates soft modes for different k-vectors; supercell size limited by maximum number of atoms maxAtoms 
% and maximum possible increase in any given direction maxIncrease (to avoid supercells like 1x1x100 if we start with 1 atom)
% supercells_sizes describes the supercells for corresp. freqs (aka 2x1x2, etc)

% added in 9.2.1: supercells are sorted by the length of the main diagonal, to check the smallest and bulkiest one first

global ORG_STRUC
atomType  = ORG_STRUC.atomType;
goodBond  = ORG_STRUC.goodBonds;
N_val     = ORG_STRUC.NvalElectrons;
  val     = ORG_STRUC.valences;
dimension = ORG_STRUC.dimension;

coords0 = coords;
lat0 = lat;
numIons0 = numIons;
N = sum(numIons);

nMax = floor(maxAtoms/N);
nMax = min(nMax, maxIncrease);

supercells = zeros(1,3);
diagonals = 0;
% build all possible supercells
for i = 1 : nMax
 for j = 1 : nMax
  for k = 1 : nMax
    if (i*j*k*N > maxAtoms) || (i*j*k>nMax) %
      continue;
    end
    supercells = vertcat(supercells, [i j k]);
    diagonals = vertcat(diagonals, i^2 + j^2 + k^2);
  end
 end
end
supercells(1,:) = []; % remove those 0
diagonals(1) = []; % remove that 0
% sort supercells by their main diagonal
[tmp, IXd] = sort(diagonals);

all_freqs = 0;
supercells_sizes = [0 0 0];
all_eigvectors = zeros(sum(3*numIons0),1);

% calculate the frequencies/vectors for supercells
%maxCells = min(ORG_STRUC.populationSize, length(diagonals))
maxCells = length(diagonals);
for ind = 1 : maxCells
    i = supercells(IXd(ind), 1);
    j = supercells(IXd(ind), 2);
    k = supercells(IXd(ind), 3);

    %  supercell iXjXk
    kVector = zeros(1,3);
    kVector1 = zeros(1,3);
    if i > 1
      kVector(1) = 1/i;
    end
    if j > 1
      kVector(2) = 1/j;
    end
    if k > 1
      kVector(3) = 1/k;
    end
    if (i > 1) && (j > 1) && (k > 1) % 4 different sign combinations
      Nminus = 4;
      signs = [1 1 1; -1 1 1; 1 -1 1; 1 1 -1];
    elseif (i > 1) && (j > 1) % 2 different sign combinations
      Nminus = 2;
      signs = [1 1 0; -1 1 0];
    elseif (j > 1) && (k > 1) % 2 different sign combinations
      Nminus = 2;
      signs = [0 1 1; 0 -1 1];
    elseif (i > 1) && (k > 1) % 2 different sign combinations
      Nminus = 2;
      signs = [1 0 1; -1 0 1];
    else
      Nminus = 1;
      signs = [1 1 1];
    end
    for m = 1 : Nminus
      kVector1 = kVector.*signs(m,:);
      [freq, eigvector] = calcSoftModes(lat0, coords0, numIons0, atomType, ...
                                   goodBond, N_val, val, dimension, kVector1);
      all_freqs      = [all_freqs; freq];              %horizontal
      all_eigvectors = [all_eigvectors, eigvector];    %vertical
      supercells_sizes = [supercells_sizes; repmat([i j k],[length(freq), 1]) ];  
    end
end

all_eigvectors(:,1) = []; % remove those 0
all_freqs(1,:) = []; % remove those 0
supercells_sizes(1,:) = []; % remove those 0

%experiment: let's sort
[all_freqs, IX]  = sort(all_freqs);
all_eigvectors   = all_eigvectors(:,IX);
supercells_sizes = supercells_sizes(IX,:);

