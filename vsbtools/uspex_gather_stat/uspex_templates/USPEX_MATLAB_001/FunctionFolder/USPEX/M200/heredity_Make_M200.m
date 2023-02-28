function [numIons, offspring, Lattice, frac, dim] = heredity_make_M200(par_one, par_two)

global POP_STRUC
global ORG_STRUC

corr        = POP_STRUC.correlation_coefficient;
cor_dir     = POP_STRUC.cor_dir;
atomType    = ORG_STRUC.atomType;
maxAt       = ORG_STRUC.maxAt;
minAt       = ORG_STRUC.minAt;
constLat    = ORG_STRUC.constLattice;

% Initiallization of parents
parent1  = POP_STRUC.POPULATION(par_one).COORDINATES;
lat1     = POP_STRUC.POPULATION(par_one).LATTICE;
order1   = POP_STRUC.POPULATION(par_one).order;
numIons1 = POP_STRUC.POPULATION(par_one).numIons;
parent2  = POP_STRUC.POPULATION(par_two).COORDINATES;
lat2     = POP_STRUC.POPULATION(par_two).LATTICE;
order2   = POP_STRUC.POPULATION(par_two).order;
numIons2 = POP_STRUC.POPULATION(par_two).numIons;

% Output: fracFrac determines what spatial fraction from parent 
frac = 0.25 + rand(1)/2;   % the most extreme difference is 1:3

% Output: choose which direct to cut
dim  = RandInt(1,1,[1,2]);     %choose only x-y directions

% Make random shift for both parents
parent1 = bsxfun(@plus, parent1, rand(1,3));
parent2 = bsxfun(@plus, parent2, rand(1,3)); 
parent1 = parent1 - floor(parent1);
parent2 = parent2 - floor(parent2);

% Choose slabs contained more/less orderer atoms
Nslab1  = heredity_Nslab(lat1, numIons1, corr, dim);
Nslab2  = heredity_Nslab(lat1, numIons2, corr, dim);
[parent1, whichIons1] = heredity_Slab(parent1, order1,...
                  Nslab1, frac, dim, cor_dir, 1);
[parent2, whichIons2] = heredity_Slab(parent2, order2,...
                  Nslab2, frac, dim, cor_dir, 2);

% Output: numIons
if ORG_STRUC.maxAt > ORG_STRUC.minAt
    N_block1 = round(numIons1/ORG_STRUC.numIons);
    N_block2 = round(numIons1/ORG_STRUC.numIons);
    numIons = RandInt(1,1,[N_block1, N_block2])*ORG_STRUC.numIons;
else
    numIons = numIons1;
end


% Output:  coordinates
offspring = heredity_coor(parent1, whichIons1, order1, numIons1,...
                          parent2, whichIons2, order2, numIons2,...
                        numIons, frac, dim, atomType, cor_dir);

% Output:  Lattice
Lattice = heredity_lattice(lat1, lat2, frac, constLat);
