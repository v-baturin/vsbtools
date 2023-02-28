function [StruNum, ParetoFront, Properties, MSG] = ParetoRankingAll()

global ORG_STRUC
global USPEX_STRUC

paretoRanking = ORG_STRUC.paretoRanking;

atomType = ORG_STRUC.atomType;
Properties = [];
%%-------------------------------------------
for i=1:length(USPEX_STRUC.POPULATION)
     gen(i)       = USPEX_STRUC.POPULATION(i).gen;
    symg(i)       = USPEX_STRUC.POPULATION(i).symg;
    enthalpy(i)   = USPEX_STRUC.POPULATION(i).Enthalpies(end);
    num(i,:)      = USPEX_STRUC.POPULATION(i).numIons;
    fit(i)        = USPEX_STRUC.POPULATION(i).Fitness;
    KPOINTS(i,:)  = USPEX_STRUC.POPULATION(i).K_POINTS(end,:);
    volume(i)     = USPEX_STRUC.POPULATION(i).Vol;
    density(i)    = USPEX_STRUC.POPULATION(i).density;
    entropy(i)    = USPEX_STRUC.POPULATION(i).struc_entr;
          s(i)    = USPEX_STRUC.POPULATION(i).S_order;
         order    = USPEX_STRUC.POPULATION(i).order;
       enth(i)    = enthalpy(i)/sum(num(i,:));
        a_o(i)    = sum(order)/sum(num(i,:));
    if ORG_STRUC.varcomp == 1
        convexhull(i)  = USPEX_STRUC.POPULATION(i).Dis_ConvexHull;
    else
        convexhull(i)  = USPEX_STRUC.POPULATION(i).Enthalpies(end)/sum(USPEX_STRUC.POPULATION(i).numIons);
    end
%    if ORG_STRUC.spin == 1
%        magmom(i,:) = sum(USPEX_STRUC.POPULATION(i).magmom_ions(end,2:end));
%        magType(i)  = USPEX_STRUC.POPULATION(i).magmom_ions(end,1);
%    else
%        magmom   =[];
%        magType  =[];
%    end
    %------------------
end
if ORG_STRUC.varcomp == 0
    MINEN  = min(convexhull);
    convexhull = convexhull - MINEN;
end
if paretoRanking ~= 1
    Properties = CalcPropforWritingPareto(paretoRanking, atomType);
end
[StruNum, ParetoFront, MSG] = PARETORANKING(fit, convexhull, Properties, paretoRanking, ORG_STRUC.optType);
if paretoRanking > 1
    Properties = -1 * Properties;
end 
%-------------QARBAL----------------------
%TOTALData = [];
%for a = 1 : length(StruNum)
%    TOTALData(a,:) = [ParetoFront(a), StruNum(a), fit(StruNum(a)), convexhull(StruNum(a))];
%end

%for a = length(StruNum): -1 : 1
%    if abs(TOTALData(a,3)) > 9999 
%	TOTALData(a,:) = [];
%    elseif TOTALData(a,4) > 1
%        TOTALData(a,:) = [];
%    elseif ORG_STRUC.optType == 6 && USPEX_STRUC.POPULATION(StruNum(a)).hardness < 2
%        TOTALData(a,:) = [];
%    end
%end
%
%for a = size(TOTALData, 1) : -1 : 2
%    for b = 1 : a - 1
%        if TOTALData(a,3) == TOTALData(b,3) && TOTALData(a,4) == TOTALData(b,4) 
%	    TOTALData(a,:) = [];
%	    break;
%        end
%   end
%end
%if length(TOTALData) ~= 0
%  ParetoFront = TOTALData(:,1);
%  StruNum     = TOTALData(:,2);
%end

%%----------------------------------------
function Prop =  CalcPropforWritingPareto(pareto, atomType)

%%% In this function we gather 3th property 
%%% for Pareto ranking.

global USPEX_STRUC

if pareto < 0
opt_sign = -1;
else
opt_sign =  1;
end
pareto = abs(pareto);
Prop = zeros(1, length(USPEX_STRUC.POPULATION));
for i = 1 : length(USPEX_STRUC.POPULATION)
    numIons = USPEX_STRUC.POPULATION(i).numIons;
    Volume = det(USPEX_STRUC.POPULATION(i).LATTICE);
    if pareto == 2 % volume
        Prop(i) = Volume/sum(numIons);
    elseif pareto == 3 % hardness
        if ~isempty(USPEX_STRUC.POPULATION(i).hardness)
            Prop(i) = -1*USPEX_STRUC.POPULATION(i).hardness;
        else
            Prop(i) = 100000;
        end
    elseif pareto == 4 % Structural Order
        if ~isempty(USPEX_STRUC.POPULATION(i).S_order)
            Prop(i) = -1*USPEX_STRUC.POPULATION(i).S_order;
        else
            Prop(i) = 100000;
        end
    elseif pareto == 5  % density
        Prop(i) = -1*calcDensity(numIons, atomType, Volume);
    elseif pareto == 6 % dielectric tensor (or rather Tr()/3)
        Prop(i) = -1*sum(USPEX_STRUC.POPULATION(i).dielectric_tensor(1:3))/3;
    elseif pareto == 7 % bandgap
        Prop(i) = -1*USPEX_STRUC.POPULATION(i).gap;
    elseif pareto == 8 % Tr(dielectric tensor)/3 multiplied by gap^2 value
        Egc = 4; %critical bandgap for semiconductor and insulator, eV
        if USPEX_STRUC.POPULATION(i).gap >= Egc
            Prop(i) = -1*(sum(USPEX_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(USPEX_STRUC.POPULATION(i).gap/Egc)^2; %Insulator
        else
            Prop(i) = -1*(sum(USPEX_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(USPEX_STRUC.POPULATION(i).gap/Egc)^6; %Semiconductor
        end
    elseif pareto == 9
        Prop(i) = -1*USPEX_STRUC.POPULATION(i).mag_moment/Volume;
    elseif pareto == 10 % structure entropy
        if ~isempty(USPEX_STRUC.POPULATION(i).struc_entr)
            Prop(i) = -1*USPEX_STRUC.POPULATION(i).struc_entr;
        else
            Prop(i) = 100000;
        end
    elseif pareto == 11
        Prop(i) = -1*USPEX_STRUC.POPULATION(i).birefringence;
    elseif pareto == 14 % thermoelectric property, replaced power factor
        Prop(i) = -1*USPEX_STRUC.POPULATION(i).TE_property;
    elseif (pareto > 1100) && (pareto < 1112)
        whichPara= mod(pareto,110);
        for i = 1 : length(USPEX_STRUC.POPULATION)
            if isempty(USPEX_STRUC.POPULATION(i).elasticProperties) || (USPEX_STRUC.POPULATION(i).elasticProperties(end)==0)
                Prop(i)=NaN;
            else
                Prop(i) = -1*USPEX_STRUC.POPULATION(i).elasticProperties(whichPara);
            end
        end
    end            
end

Prop = opt_sign*Prop; % change the mode of optimization (minimization <=> maximization)

%%%---------------------------------------------------------------------
function [Ranking, ParetoFront, MSG] = PARETORANKING(fitness, Enthalpies, Properties, paretoRanking, optType)
global USPEX_STRUC

LENFIT = length(fitness);
MSG = 0 ;
%% first we rank properties with respect to one property and then we take the structures which 
%% pass limitations and rank them according to two properties, and then put the structures 
%% which couldn't pass the the limitations, at the end of this structures in ranking.

%% First normal ranking according to fitness
[nothing, firstRank]  = sort(fitness);
firstRank = firstRank';
Ranking   = firstRank;
ParetoF1(:,1) = 1 : LENFIT;
ParetoFront = ParetoF1;
firstStep = [Ranking, ParetoFront];

ii = 0;
for i = 1 : LENFIT
    if Enthalpies(i) < 1 && abs(fitness(i)) < 9999  
        ii = ii + 1 ;
        ID(ii)   = i;
        ENE(ii)  = Enthalpies(i);
        Fit(ii)  = fitness(i);
        if paretoRanking ~= 1
  	    PROP(ii) = Properties(i);
        end	
    end
end
%% Second Rank according to 2 or 3 properties
if paretoRanking == 1
    [secondRank, Dist, ParetoF2] = ParetoRanking_2D(Fit, ENE, ID);    % IDs are the number of Good Structures
    secondStep = [secondRank, ParetoF2];
    [Ranking, ParetoFront] = RANKING(secondStep, firstStep);
else
    [secondRank, Dist, ParetoF2] = ParetoRanking_3D(Fit, ENE, PROP, ID);
    secondStep = [secondRank, ParetoF2];
    [Ranking, ParetoFront] = RANKING(secondStep, firstStep);
end
initialStep = [Ranking, ParetoFront];
%% Third: Only for Dielectric constant, again ranking 
%% with respect to 2 or 3 properties, this time we have one 
%% more constrain, if no structures satisfy these constrains
%% ranking will not be performed, and previous rank is used.
if optType == 6
    Fit = []; ENE = []; ID = []; PROP = [];
    ii      = 0;
    for i = 1 : LENFIT
       if Enthalpies(i) < 1 && abs(fitness(i)) < 9999 && fitness(i) > 1 && USPEX_STRUC.POPULATION(i).hardness > 2  %% Here we don't consider stability
            ii = ii + 1 ;
            ID(ii)   = i;
            ENE(ii)  = Enthalpies(i);
            Fit(ii)  = fitness(i);
            if paretoRanking ~= 1
                PROP(ii) = Properties(i);
            end
        end
    end
    if ii ~= 0
        if paretoRanking == 1
	    [thirdRank, Dist, ParetoF3] = ParetoRanking_2D(Fit, ENE, ID);
            thirdStep = [thirdRank, ParetoF3];
            [Ranking, ParetoFront] = RANKING(thirdStep, initialStep);
        else
            [thirdRank, Dist, ParetoF3] = ParetoRanking_3D(Fit, ENE, PROP, ID);
            thirdStep = [thirdRank, ParetoF3];
            [Ranking, ParetoFront] = RANKING(thirdStep, initialStep);
        end
    else
	disp('This system cannot satisfy hardness > 2GPa condition, ranking is performed without this constrain!');
	MSG   = 1 ;
    end
end

%%%----------------------------------------------
function [Ranking, ParetoFront] = RANKING(BestRank, WorseRank)

if length(BestRank) ~= 0
   for a = size(WorseRank,1) : -1 : 1
     for b = 1 : size(BestRank,1)
        if WorseRank(a,1) == BestRank(b,1)
            WorseRank(a,:) = [];
            break;
        end
     end
   end
   if length(WorseRank) ~= 0
      WorseRank(:,2) = WorseRank(:,2) + BestRank(end,2);
      if WorseRank(1,2) ~= (BestRank(end,2) + 1)
         Diff = WorseRank(1,2) - (BestRank(end,2) + 1);
         WorseRank(:,2) = WorseRank(:,2) - Diff;
      end
      for a = 2 : size(WorseRank,1)
         if WorseRank(a,2) > WorseRank(a-1,2) + 1
	    Diff = WorseRank(a,2) - (WorseRank(a-1,2) + 1);
            WorseRank(a:end,2) = WorseRank(a:end,2) - Diff;
         end
      end
   end
   result = [BestRank ; WorseRank];
   Ranking     = result(:,1);
   ParetoFront = result(:,2);
else
   Ranking     = WorseRank(:,1);
   ParetoFront = WorseRank(:,2);
end

