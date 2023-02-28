function MINHOP_SoftModeMutation(Ind_No, ind, the_same_structure)

global POP_STRUC
global ORG_STRUC
global OFF_STRUC

atomType  = ORG_STRUC.atomType;
dimension = ORG_STRUC.dimension;
minDist   = ORG_STRUC.minDistMatrice;
goodBond  = ORG_STRUC.goodBonds;
N_val     = ORG_STRUC.NvalElectrons;
  val     = ORG_STRUC.valences;


mutationRate  =  ORG_STRUC.mutationRate;
goodAtomMutant = 0;
structureFailed = 0;
goodMutLattice = 0;
superCellSize = [1 1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%be carefull%     lat = OFF_STRUC.lat;


  N = sum(POP_STRUC.POPULATION(ind).numIons);
  lat  = POP_STRUC.POPULATION(ind).LATTICE;
  coor = POP_STRUC.POPULATION(ind).COORDINATES;
  numIons = POP_STRUC.POPULATION(ind).numIons;
  [freq, eigvector] = calcSoftModes(lat, coor, numIons, atomType, goodBond, N_val, val, dimension);

while goodAtomMutant ~= 1
       logic1=0;
%%%%%%%%%%Linear combination of modes (by defualt 10 percent of modes)
 natom = sum(numIons);
 mode = 4;
 lcom = round(0.3*natom) + mode;
 normali_f = 0;
 stoch_vect = zeros(3*natom,1);
 superCellSize = [1 1 1];
 while  mode <= lcom  % stochastic softmutation
    eigV = eigvector(:,mode);
    if isAcoustic(eigV, coor, superCellSize)
            mode = mode + 1;
            lcom = lcom +1;
            if natom > 1
              continue
            elseif sum(superCellSize) == 1 % 1 atom, translational accoustic mode
              continue
            end
    else
    c = 2*rand-1;
    stoch_vect = stoch_vect + c*eigvector(:,mode);
    mode = mode + 1;
    normali_f = normali_f + c^2;
    end
 end

%%%%%%%%%%             

             coord0 = coor; 
             numIons0 = numIons;
             lat0 = lat;
             superCellSize = [1 1 1];

if POP_STRUC.generation < 2 || the_same_structure == 0
    last_mut = 1;
else
last_mut = POP_STRUC.POPULATION(ind).last_mut;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Start Minima Hopping  %%%%%%%%%%%%%%%%%%%%%%%%%
          i = last_mut;
           [MUT_LAT, MUT_COORD, deviation] = move_along_SoftMode_minima_hopping(coord0, numIons0, lat0, stoch_vect/sqrt(normali_f), 1+i/20, (i+1)*mutationRate);
%%%%%To Optimize the lattice
if ORG_STRUC.constLattice ~= 1
    coordmah = MUT_COORD*MUT_LAT;
    [coordmah, MUT_LAT] = optLattice(coordmah, MUT_LAT);
    MUT_COORD = coordmah/MUT_LAT;
end
%%%%%To optimize the lattice
           goodAtomMutant = distanceCheck(MUT_COORD, MUT_LAT, numIons0, ORG_STRUC.minDistMatrice);
           structureFailed = structureFailed + 1 - goodAtomMutant ;
           if structureFailed >= 100
             disp('Distance check failed more than 100 times in a single softmutation.')
             structureFailed = 0;
             logic1 = 0;
           end
           if goodAtomMutant 

              logic1 = 1;  
              if the_same_structure
                  last_mut = last_mut + 1;
              end
              OFF_STRUC.POPULATION(ind).COORDINATES = MUT_COORD;
              OFF_STRUC.POPULATION(ind).LATTICE = MUT_LAT;
              OFF_STRUC.lat = MUT_LAT;
              POP_STRUC.POPULATION(ind).COORDINATES = MUT_COORD;
              POP_STRUC.POPULATION(ind).LATTICE = MUT_LAT;
              POP_STRUC.POPULATION(ind).superCell = superCellSize.*POP_STRUC.POPULATION(ind).superCell;
              info_parents = struct('parent', {},'mut_degree', {},'mut_mode',{},'mut_fre',{},'enthalpy',{});
              info_parents(1).parent = ind;
              info_parents.mut_degree = deviation;
              info_parents.enthalpy = POP_STRUC.POPULATION(ind).Enthalpies(end);
        disp(['Structure softmutated with combination of modes for the # ' num2str(last_mut) 'th time); Mut_degree_max = ' num2str(max(deviation)) ', Mut_degree_aver = ' num2str(mean(deviation))]);

              POP_STRUC.POPULATION(ind).Parents = info_parents;
              POP_STRUC.POPULATION(ind).numIons = numIons0;
              OFF_STRUC.POPULATION(ind).last_mut = last_mut;

              break;
           elseif (i == 100)
              break;         
           end
 
end
 if (logic1 == 0) 
   % atomMutation(ind);
    disp('We should figure out what to do if method fails after certain iteration for example 100 cycle');
    goodAtomMutant = 1;
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%   END Minima Hopping   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
