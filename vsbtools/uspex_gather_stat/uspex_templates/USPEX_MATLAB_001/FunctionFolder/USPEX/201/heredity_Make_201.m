function [numIons, offspring, Lattice, frac, dim] = heredity_make_201(par_one, par_two, cell)

global POP_STRUC
global ORG_STRUC

corr        = POP_STRUC.correlation_coefficient;
cor_dir     = POP_STRUC.cor_dir;
atomType    = ORG_STRUC.atomType;
lat_base    = ORG_STRUC.bulk_lat;

% Initiallization of parents
parent1  = POP_STRUC.POPULATION(par_one).COORDINATES;
lat1     = POP_STRUC.POPULATION(par_one).LATTICE;
order1   = POP_STRUC.POPULATION(par_one).order;
numIons1 = POP_STRUC.POPULATION(par_one).numIons;
parent2  = POP_STRUC.POPULATION(par_two).COORDINATES;
lat2     = POP_STRUC.POPULATION(par_two).LATTICE;
order2   = POP_STRUC.POPULATION(par_two).order;
numIons2 = POP_STRUC.POPULATION(par_two).numIons;
[lat1, parent1, order1, numIons1]=cresupercell(lat1, parent1, numIons1, lat_base, cell, atomType);
[lat2, parent2, order2, numIons2]=cresupercell(lat2, parent2, numIons2, lat_base, cell, atomType);
%disp('heredity_make_201')

% Output: fracFrac determines what spatial fraction from parent 
frac = 0.25 + rand(1)/2;   % the most extreme difference is 1:3

% Output: choose which direct to cut
dim  = RandInt(1,1,[1,2]);     %choose only x-y directions

% Make random shift for both parents
parent1 = bsxfun(@plus, parent1, [rand(1,2) 0]);
parent2 = bsxfun(@plus, parent2, [rand(1,2) 0]); 
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
numIons = round((numIons1+numIons2)/2);


% Output:  coordinates
offspring = heredity_coor(parent1, whichIons1, order1, numIons1,...
                          parent2, whichIons2, order2, numIons2,...
                        numIons, frac, dim, atomType, cor_dir);

% Output:  Lattice
Lattice = lat1;
