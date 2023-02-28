function fitness = CalcFitness_001()

global POP_STRUC
global ORG_STRUC

% calculation CLUSTERS structure
global CLUSTERS
load CLUSTERS.mat

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------%
%------------------- case of two-component (binary) clusters ---------------------%
%---------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(ORG_STRUC.numIons,2) == 2

minions1 = ORG_STRUC.numIons(1,1);
minions2 = ORG_STRUC.numIons(1,2);

for i = 1:length(POP_STRUC.POPULATION)
    cx = POP_STRUC.POPULATION(i).numIons(1);
    cy = POP_STRUC.POPULATION(i).numIons(2);
    cxx = cx-minions1+2;
    cyy = cy-minions2+2;
    if POP_STRUC.POPULATION(i).Enthalpies(end) <= CLUSTERS.composition(cxx,cyy).bestEnthalpy
        CLUSTERS.composition(cxx,cyy).bestEnthalpy = POP_STRUC.POPULATION(i).Enthalpies(end);
	CLUSTERS.composition(cxx,cyy).origin = i;
	CLUSTERS.composition(cxx,cyy).COORDINATES = POP_STRUC.POPULATION(i).COORDINATES;
	CLUSTERS.composition(cxx,cyy).numIons     = POP_STRUC.POPULATION(i).numIons;
	CLUSTERS.composition(cxx,cyy).LATTICE     = POP_STRUC.POPULATION(i).LATTICE;
	CLUSTERS.composition(cxx,cyy).symg        = POP_STRUC.POPULATION(i).symg;
	CLUSTERS.composition(cxx,cyy).order       = POP_STRUC.POPULATION(i).order;
	CLUSTERS.composition(cxx,cyy).Number      = POP_STRUC.POPULATION(i).Number;
    end    
end

% calculation CLUSTERS.ConvexHall
numcomp = 0;
CLUSTERS.ConvexHall = [];
for i = ORG_STRUC.numIons(1,1) : ORG_STRUC.numIons(2,1)
for j = ORG_STRUC.numIons(1,2) : ORG_STRUC.numIons(2,2)
   ii = i-minions1+2;
   jj = j-minions2+2;
   d2Ex = CLUSTERS.composition(ii+1,jj).bestEnthalpy + CLUSTERS.composition(ii-1,jj).bestEnthalpy - 2*CLUSTERS.composition(ii,jj).bestEnthalpy;
   d2Ey = CLUSTERS.composition(ii,jj+1).bestEnthalpy + CLUSTERS.composition(ii,jj-1).bestEnthalpy - 2*CLUSTERS.composition(ii,jj).bestEnthalpy;
   if (d2Ex >= 0) && (d2Ey >= 0)
      numcomp = numcomp + 1;
      CLUSTERS.ConvexHall(numcomp).Enthalpy = CLUSTERS.composition(ii,jj).bestEnthalpy;
      CLUSTERS.ConvexHall(numcomp).numIons = [i, j];
      CLUSTERS.ConvexHall(numcomp).origin = CLUSTERS.composition(ii,jj).origin;
   end
end
end
CLUSTERS.number_ConvexHall = numcomp;
% adding one more fictitious field with compostion (1000,1001) to avoid error 
% when considering binary clusters when for one element the number of atoms is fixed
if (ORG_STRUC.numIons(2,1)==ORG_STRUC.numIons(1,1)) || (ORG_STRUC.numIons(2,2)==ORG_STRUC.numIons(1,2))
   CLUSTERS.ConvexHall(end+1).Enthalpy = 0;
   CLUSTERS.ConvexHall(end).numIons = [1000, 1001];
end

%safesave ('CLUSTERS.mat', CLUSTERS)

% calculation POP_STRUC.convex_hull
POP_STRUC.convex_hull = [];
for i = 1 : numcomp
   POP_STRUC.convex_hull(i,1) = CLUSTERS.ConvexHall(i).numIons(1);
   POP_STRUC.convex_hull(i,2) = CLUSTERS.ConvexHall(i).numIons(2);
   POP_STRUC.convex_hull(i,3) = CLUSTERS.ConvexHall(i).Enthalpy/(CLUSTERS.ConvexHall(i).numIons(1)+CLUSTERS.ConvexHall(i).numIons(2));
   POP_STRUC.convex_hull(i,4) = CLUSTERS.ConvexHall(i).origin;
end

% calculation fitness
fitness = zeros(1,length(POP_STRUC.POPULATION));

calc_type = 'mag';

if (calc_type == 'all')
    for fit_loop = 1:length(POP_STRUC.POPULATION)
        cx = POP_STRUC.POPULATION(fit_loop).numIons(1);
        cy = POP_STRUC.POPULATION(fit_loop).numIons(2);
        cxx = cx-minions1+2;
        cyy = cy-minions2+2;
        fitness(fit_loop) = (POP_STRUC.POPULATION(fit_loop).Enthalpies(end) - CLUSTERS.composition(cxx,cyy).bestEnthalpy);
    end
elseif (calc_type == 'mag')
    for fit_loop = 1:length(POP_STRUC.POPULATION)
        cx = POP_STRUC.POPULATION(fit_loop).numIons(1);
        cy = POP_STRUC.POPULATION(fit_loop).numIons(2);
	done = 0;
	for i = 1:numcomp
	    x_ch = CLUSTERS.ConvexHall(i).numIons(1);
	    y_ch = CLUSTERS.ConvexHall(i).numIons(2);
	    if (cx==x_ch) && (cy==y_ch)
		done = 1;
		Eref = CLUSTERS.ConvexHall(i).Enthalpy;
	    end
	end
	
	if ~done
	    for i=1:length(CLUSTERS.ConvexHall) %numcomp
		xch = CLUSTERS.ConvexHall(i).numIons(1);
		ych = CLUSTERS.ConvexHall(i).numIons(2);
		r(i) = sqrt((cx-xch)^2 + (cy-ych)^2);
	    end
       	    [r_upor, n_upor] = sort(r);
	    n_current = 2;
	    while (~done)
		n_current = n_current + 1;
		for i = 1:n_current-2
		    for j = i+1:n_current-1
			x1 = CLUSTERS.ConvexHall(n_upor(i)).numIons(1);
	    		y1 = CLUSTERS.ConvexHall(n_upor(i)).numIons(2);
	    		x2 = CLUSTERS.ConvexHall(n_upor(j)).numIons(1);
	    		y2 = CLUSTERS.ConvexHall(n_upor(j)).numIons(2);
			x3 = CLUSTERS.ConvexHall(n_upor(n_current)).numIons(1);
	    		y3 = CLUSTERS.ConvexHall(n_upor(n_current)).numIons(2);
	    		S = AreaTriangle(x1,y1,x2,y2,x3,y3);
	    		S1 = AreaTriangle(cx,cy,x1,y1,x2,y2);
	    		S2 = AreaTriangle(cx,cy,x1,y1,x3,y3);
	    		S3 = AreaTriangle(cx,cy,x2,y2,x3,y3);
			if (S~=0) && (S1+S2+S3 <= S+0.00000001)
			    E1 = CLUSTERS.ConvexHall(n_upor(i)).Enthalpy;
			    E2 = CLUSTERS.ConvexHall(n_upor(j)).Enthalpy;
			    E3 = CLUSTERS.ConvexHall(n_upor(n_current)).Enthalpy;
			    Eref=(det([x1 y1 E1;x2 y2 E2;x3 y3 E3])-cx*det([y1 E1 1;y2 E2 1;y3 E3 1])+cy*det([x1 E1 1;x2 E2 1;x3 E3 1]))/det([x1 y1 1;x2 y2 1;x3 y3 1]);
			    done = 1;
			end
		    if done break; end
		    end
		if done break; end
		end
	    end

	end
    fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end) - Eref;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------%
%----------------------- case of three-component clusters --------------------------%
%---------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(ORG_STRUC.numIons,2) == 3

minions1 = ORG_STRUC.numIons(1,1);
minions2 = ORG_STRUC.numIons(1,2);
minions3 = ORG_STRUC.numIons(1,3);

for i = 1:length(POP_STRUC.POPULATION)
    cx = POP_STRUC.POPULATION(i).numIons(1);
    cy = POP_STRUC.POPULATION(i).numIons(2);
    cz = POP_STRUC.POPULATION(i).numIons(3);
    cxx = cx-minions1+2;
    cyy = cy-minions2+2;
    czz = cz-minions3+2;
    
    if POP_STRUC.POPULATION(i).Enthalpies(end) <= CLUSTERS.composition(cxx,cyy,czz).bestEnthalpy
        CLUSTERS.composition(cxx,cyy,czz).bestEnthalpy = POP_STRUC.POPULATION(i).Enthalpies(end);
        CLUSTERS.composition(cxx,cyy,czz).origin = i;
        CLUSTERS.composition(cxx,cyy,czz).COORDINATES = POP_STRUC.POPULATION(i).COORDINATES;
        CLUSTERS.composition(cxx,cyy,czz).numIons     = POP_STRUC.POPULATION(i).numIons;
        CLUSTERS.composition(cxx,cyy,czz).LATTICE     = POP_STRUC.POPULATION(i).LATTICE;
        CLUSTERS.composition(cxx,cyy,czz).symg        = POP_STRUC.POPULATION(i).symg;
        CLUSTERS.composition(cxx,cyy,czz).order       = POP_STRUC.POPULATION(i).order;
        CLUSTERS.composition(cxx,cyy,czz).Number      = POP_STRUC.POPULATION(i).Number;
    end    
end

% calculation CLUSTERS.ConvexHall
numcomp = 0;
CLUSTERS.ConvexHall = [];
for i = ORG_STRUC.numIons(1,1) : ORG_STRUC.numIons(2,1)
for j = ORG_STRUC.numIons(1,2) : ORG_STRUC.numIons(2,2)
for k = ORG_STRUC.numIons(1,3) : ORG_STRUC.numIons(2,3)
   ii = i-minions1+2;
   jj = j-minions2+2;
   kk = k-minions3+2;
   
   d2Ex = CLUSTERS.composition(ii+1,jj,kk).bestEnthalpy + CLUSTERS.composition(ii-1,jj,kk).bestEnthalpy - 2*CLUSTERS.composition(ii,jj,kk).bestEnthalpy;
   d2Ey = CLUSTERS.composition(ii,jj+1,kk).bestEnthalpy + CLUSTERS.composition(ii,jj-1,kk).bestEnthalpy - 2*CLUSTERS.composition(ii,jj,kk).bestEnthalpy;
   d2Ez = CLUSTERS.composition(ii,jj,kk+1).bestEnthalpy + CLUSTERS.composition(ii,jj,kk-1).bestEnthalpy - 2*CLUSTERS.composition(ii,jj,kk).bestEnthalpy;
   
   if (d2Ex >= 0) && (d2Ey >= 0) && (d2Ez >= 0)
      numcomp = numcomp + 1;
      CLUSTERS.ConvexHall(numcomp).Enthalpy = CLUSTERS.composition(ii,jj,kk).bestEnthalpy;
      CLUSTERS.ConvexHall(numcomp).numIons = [i,j,k];
      CLUSTERS.ConvexHall(numcomp).origin = CLUSTERS.composition(ii,jj,kk).origin;
   end
end
end
end
CLUSTERS.number_ConvexHall = numcomp;

% % adding one more fictitious field with compostion (1000,1001) to avoid error 
% % when considering binary clusters when for one element the number of atoms is fixed
% if (ORG_STRUC.numIons(2,1)==ORG_STRUC.numIons(1,1)) || (ORG_STRUC.numIons(2,2)==ORG_STRUC.numIons(1,2))
%    CLUSTERS.ConvexHall(end+1).Enthalpy = 0;
%    CLUSTERS.ConvexHall(end).numIons = [1000, 1001];
% end

%safesave ('CLUSTERS.mat', CLUSTERS)

% calculation POP_STRUC.convex_hull
POP_STRUC.convex_hull = [];
for i = 1 : numcomp
   POP_STRUC.convex_hull(i,1) = CLUSTERS.ConvexHall(i).numIons(1);
   POP_STRUC.convex_hull(i,2) = CLUSTERS.ConvexHall(i).numIons(2);
   POP_STRUC.convex_hull(i,3) = CLUSTERS.ConvexHall(i).numIons(3);
   POP_STRUC.convex_hull(i,4) = CLUSTERS.ConvexHall(i).Enthalpy/sum(CLUSTERS.ConvexHall(i).numIons);
   POP_STRUC.convex_hull(i,5) = CLUSTERS.ConvexHall(i).origin;
end

% calculation fitness
fitness = zeros(1,length(POP_STRUC.POPULATION));

for fit_loop = 1:length(POP_STRUC.POPULATION)
    cx = POP_STRUC.POPULATION(fit_loop).numIons(1);
    cy = POP_STRUC.POPULATION(fit_loop).numIons(2);
    cz = POP_STRUC.POPULATION(fit_loop).numIons(3);
    cxx = cx-minions1+2;
    cyy = cy-minions2+2;
    czz = cz-minions3+2;
    
%     %all
%     fitness(fit_loop) = (POP_STRUC.POPULATION(fit_loop).Enthalpies(end) - CLUSTERS.composition(cxx,cyy,czz).bestEnthalpy);
    
    %mag
    Ex_mean = 0.5*(CLUSTERS.composition(cxx+1,cyy,czz).bestEnthalpy + CLUSTERS.composition(cxx-1,cyy,czz).bestEnthalpy);
    Ey_mean = 0.5*(CLUSTERS.composition(cxx,cyy+1,czz).bestEnthalpy + CLUSTERS.composition(cxx,cyy-1,czz).bestEnthalpy);
    Ez_mean = 0.5*(CLUSTERS.composition(cxx,cyy,czz+1).bestEnthalpy + CLUSTERS.composition(cxx,cyy,czz-1).bestEnthalpy);
    
    Eref = min([CLUSTERS.composition(cxx,cyy,czz).bestEnthalpy, Ex_mean, Ey_mean, Ez_mean]);
    fitness(fit_loop) = POP_STRUC.POPULATION(fit_loop).Enthalpies(end) - Eref;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------------%
%----------------------- case of one-component clusters --------------------------%
%---------------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif size(ORG_STRUC.numIons,2) == 1

for i = 1:length(POP_STRUC.POPULATION)
    Comp = POP_STRUC.POPULATION(i).numIons - ORG_STRUC.numIons(1,1) + 1;
    if POP_STRUC.POPULATION(i).Enthalpies(end) <= CLUSTERS.composition(Comp).bestEnthalpy
        CLUSTERS.composition(Comp).bestEnthalpy = POP_STRUC.POPULATION(i).Enthalpies(end);
        CLUSTERS.composition(Comp).origin = i;
        CLUSTERS.composition(Comp).COORDINATES = POP_STRUC.POPULATION(i).COORDINATES;
        CLUSTERS.composition(Comp).numIons     = POP_STRUC.POPULATION(i).numIons;
        CLUSTERS.composition(Comp).LATTICE     = POP_STRUC.POPULATION(i).LATTICE;
        CLUSTERS.composition(Comp).symg        = POP_STRUC.POPULATION(i).symg;
        CLUSTERS.composition(Comp).order       = POP_STRUC.POPULATION(i).order;
        CLUSTERS.composition(Comp).Number      = POP_STRUC.POPULATION(i).Number;
    end    
end

numComp_CH = 1;
numComp = CLUSTERS.number_compositions;
CLUSTERS.ConvexHall(1).Enthalpy = CLUSTERS.composition(1).bestEnthalpy;
CLUSTERS.ConvexHall(1).numIons = CLUSTERS.composition(1).numIons;
CLUSTERS.ConvexHall(1).origin = CLUSTERS.composition(1).origin;

for i = 2:numComp-1
    if (2*CLUSTERS.composition(i).bestEnthalpy <= CLUSTERS.composition(i-1).bestEnthalpy + CLUSTERS.composition(i+1).bestEnthalpy)
        numComp_CH = numComp_CH + 1;
        CLUSTERS.ConvexHall(numComp_CH).Enthalpy = CLUSTERS.composition(i).bestEnthalpy;
        CLUSTERS.ConvexHall(numComp_CH).numIons = CLUSTERS.composition(i).numIons;
	CLUSTERS.ConvexHall(numComp_CH).origin = CLUSTERS.composition(i).origin;
    end
end

CLUSTERS.ConvexHall_numComp = numComp_CH + 1;
CLUSTERS.ConvexHall(numComp_CH + 1).Enthalpy = CLUSTERS.composition(numComp).bestEnthalpy;
CLUSTERS.ConvexHall(numComp_CH + 1).numIons = CLUSTERS.composition(numComp).numIons;
CLUSTERS.ConvexHall(numComp_CH + 1).origin = CLUSTERS.composition(numComp).origin;

CLUSTERS.ConvexHall(numComp_CH + 2).numIons = CLUSTERS.ConvexHall(numComp_CH + 1).numIons + 1;
CLUSTERS.ConvexHall(numComp_CH + 2).Enthalpy = 0;
%safesave ('CLUSTERS.mat', CLUSTERS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fitness = zeros(1,length(POP_STRUC.POPULATION));
%Enthalpy
for fit_loop = 1:length(POP_STRUC.POPULATION)
    i=1;
    curIons = POP_STRUC.POPULATION(fit_loop).numIons;
    while 1
       i=i+1;
       if (CLUSTERS.ConvexHall(i).numIons > curIons)
          break;
       end
    end
    n1 = CLUSTERS.ConvexHall(i-1).numIons;
    n2 = CLUSTERS.ConvexHall(i).numIons;
    E1 = CLUSTERS.ConvexHall(i-1).Enthalpy;
    E2 = CLUSTERS.ConvexHall(i).Enthalpy;
    fitness(fit_loop) = (POP_STRUC.POPULATION(fit_loop).Enthalpies(end) - (E1+(curIons-n1)*(E2-E1)/(n2-n1))); % / curIons;
end

POP_STRUC.convex_hull = [];
for i = 1 : CLUSTERS.ConvexHall_numComp
   POP_STRUC.convex_hull(i,1) = CLUSTERS.ConvexHall(i).numIons;
   POP_STRUC.convex_hull(i,2) = CLUSTERS.ConvexHall(i).Enthalpy/CLUSTERS.ConvexHall(i).numIons;
   POP_STRUC.convex_hull(i,3) = CLUSTERS.ConvexHall(i).origin;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : length(fitness)
    if fitness(i)==0
        fitness(i) = -0.0001*rand(1);
    end
end

fitness = ORG_STRUC.opt_sign*fitness; % change the mode of optimization (minimization <=> maximization)
fitnessRange = FitnessLimits(ORG_STRUC.opt_sign * ORG_STRUC.optType);
for i = 1 : length(fitness)
    if POP_STRUC.POPULATION(i).Enthalpies(end) > 99999    % structure with errors
        fitness(i) = 100000;
    elseif fitness(i) < fitnessRange(1)
        fitness(i) = 100000;
    elseif fitness(i) > fitnessRange(2)
        fitness(i) = 100000;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if POP_STRUC.generation==1 
  F = [];
  for i=1:length(POP_STRUC.POPULATION)
    F(i) = fitness(i);
  end
  num_F = find(abs(F-mean(F))<std(F)); % number of fitneses which haven't border values (causes a bug in fitness graph)
  F_ots=[];
  for i=1:length(num_F)
    F_ots(i)=F(num_F(i));
  end
  CLUSTERS.fitness_range = 2*std(F_ots); %this field is diapason of y-axis in graph of fitness
end
safesave ('CLUSTERS.mat', CLUSTERS)
