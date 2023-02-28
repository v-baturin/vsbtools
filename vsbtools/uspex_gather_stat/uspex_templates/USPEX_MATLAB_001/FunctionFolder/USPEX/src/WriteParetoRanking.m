function WriteParetoRanking(resFolder)

global ORG_STRUC
global USPEX_STRUC

paretoRanking = ORG_STRUC.paretoRanking;
FITNESS = PropertiesName(ORG_STRUC.optType);
ProName = PropertiesName(paretoRanking);

%Properties = [];
fpath = [resFolder '/Pareto_ranking'];
fp = fopen(fpath, 'w+');
%%-------------------------------------------
for i=1:length(USPEX_STRUC.POPULATION)
     gen(i)       = USPEX_STRUC.POPULATION(i).gen;
    symg(i)       = USPEX_STRUC.POPULATION(i).symg;
    enthalpy(i)   = USPEX_STRUC.POPULATION(i).Enthalpies(end);
    num(i,:)      = USPEX_STRUC.POPULATION(i).numIons;
    if ORG_STRUC.optType == 2
	fit(i)    = USPEX_STRUC.POPULATION(i).Fitness;
    else
        fit(i)    = -1*ORG_STRUC.opt_sign*USPEX_STRUC.POPULATION(i).Fitness;
    end
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
%if paretoRanking ~= 1
%    Properties = CalcPropforWritingPareto(paretoRanking);
%end
%[StruNum, ParetoFront, MSG] = PARETORANKING(fit, convexhull, Properties, paretoRanking, ORG_STRUC.optType);
%if paretoRanking > 1
%    Properties = -1 * Properties;
%end 
%ID = StruNum;
%ParetoF = ParetoFront;
[StruNum, ParetoFront, Properties, MSG] = ParetoRankingAll();

ranking = StruNum;
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
for i = 2 : N
    if USPEX_STRUC.POPULATION(ranking(i)).ToCount > 0
        same = 0;
        for j= 1 : i-1
            if USPEX_STRUC.POPULATION(ranking(j)).ToCount > 0
                if SameStructure_order(ranking(i), ranking(j), USPEX_STRUC)
                    USPEX_STRUC.POPULATION(ranking(i)).ToCount = 0;
                    break;
                end
            end
        end
    end
end
%--------------------------------------------------------------------------

%-------------QARBAL----------------------
TOTALData = [];
for a = 1 : length(StruNum)
    TOTALData(a,:) = [ParetoFront(a), StruNum(a), fit(StruNum(a)), convexhull(StruNum(a))];
end

for a = length(StruNum): -1 : 1
    if abs(TOTALData(a,3)) > 9999 
	TOTALData(a,:) = [];
    elseif length(num2str(TOTALData(a,3))) == 3
	if num2str(TOTALData(a,3)) == 'NaN'
            TOTALData(a,:) = [];
	end
    elseif USPEX_STRUC.POPULATION(StruNum(a)).ToCount == 0
	TOTALData(a,:) = [];
    elseif TOTALData(a,4) > 1
        TOTALData(a,:) = [];
    elseif ORG_STRUC.optType == 6 && USPEX_STRUC.POPULATION(StruNum(a)).hardness < 2
        TOTALData(a,:) = []; 
    elseif ORG_STRUC.optType == 6 && TOTALData(a,3) < 1
        TOTALData(a,:) = [];
    elseif ORG_STRUC.optType == 3 &&  TOTALData(a,3) < 0.5
        TOTALData(a,:) = [];
    end
end

for a = size(TOTALData, 1) : -1 : 2
    for b = 1 : a - 1
        if TOTALData(a,3) == TOTALData(b,3) && TOTALData(a,4) == TOTALData(b,4) 
	    TOTALData(a,:) = [];
	    break;
        end
   end
end
if length(TOTALData) ~= 0
  ParetoFront = TOTALData(:,1);
  StruNum     = TOTALData(:,2);
  makeParetoFigures(TOTALData(:,3), TOTALData(:,4), resFolder);
end
%%----------------------------------------
if MSG == 1 
    fprintf(fp, '*** This system couldnt satisfy hardness > 2GPa condition, ranking is performed without this constrain!!! \n');
end
if paretoRanking == 1
    fprintf(fp, [ 'Pareto  ID   Origin     Composition     Enthalpy   Volume  Density ' FITNESS ' ConvexHull  KPOINTS   SYMM  Q_entr A_order S_order \n']);
else
    fprintf(fp, [ 'Pareto  ID   Origin     Composition     Enthalpy   Volume  Density ' FITNESS ' ConvexHull ' ProName ' KPOINTS   SYMM  Q_entr A_order S_order \n']);
end
fprintf(fp, ['front                                   eV/atom    (A^3)  (g/cm^3)\n']);
%%-----------------------------------------
comp_format = '[%11s]';
volume_format = '%9.3f';
density_format = '%7.3f';
fitness_format = '%12.3f';
convex_format  = '%12.3f';
property_format= '%12.3f';
kpoints_format = '[%2d %2d %2d]';
symmetry_format = '%3d';
entropy_format = '%8.3f';
ao_format = '%6.3f';
s_format = '%6.3f';
%mag_format = '%8.3f';
%magType_format = '%6s';
%%------------------------------------------
for i = 1 : length(StruNum)
    composition = sprintf('%3d',num(StruNum(i),:));
    shift=[4, 2, 1]; %so far we only consider 6 component
    if size(composition,2)<11
       composition=[composition,blanks(shift(length(num(StruNum(i),:))))];
    end
   if paretoRanking == 1
     fprintf(fp,['%3d %6d  %-11s '  comp_format ' %9.3f ' volume_format ' ' density_format ' ' fitness_format ' ' convex_format '    ' kpoints_format ' ' symmetry_format ' ' entropy_format ' ' ao_format ' ' s_format  ' \n'], ...
     ParetoFront(i) , StruNum(i), USPEX_STRUC.POPULATION(StruNum(i)).howCome, composition,  enth(StruNum(i)),   volume(StruNum(i)), density(StruNum(i)),...
     fit(StruNum(i)), convexhull(StruNum(i)), KPOINTS((StruNum(i)),:), symg(StruNum(i)), entropy(StruNum(i)), a_o(StruNum(i)), s(StruNum(i)));
   else
     fprintf(fp,['%3d %6d  %-11s '  comp_format ' %9.3f ' volume_format ' ' density_format ' ' fitness_format ' ' convex_format ' ' property_format '    ' kpoints_format ' ' symmetry_format ' ' entropy_format ' ' ao_format ' ' s_format  '\n'], ...
     ParetoFront(i) , StruNum(i), USPEX_STRUC.POPULATION(StruNum(i)).howCome, composition,  enth(StruNum(i)),   volume(StruNum(i)), density(StruNum(i)),...
     fit(StruNum(i)), convexhull(StruNum(i)), Properties(StruNum(i)),  KPOINTS((StruNum(i)),:), symg(StruNum(i)), entropy(StruNum(i)),...
     a_o(StruNum(i)), s(StruNum(i)));
   end
end
fclose(fp);
%%%%%%------------------------------------------------------------------------
function ProName = PropertiesName(OptType)

ProName = '';

if OptType == 2
    ProName = '     Volume    ';
elseif  OptType == 3
    ProName = '   Hardness    ';
elseif  OptType == 4
    ProName = '  Struc_order  ';
elseif  OptType == 5
    ProName = '    Density    ';
elseif  OptType == 6
    ProName = ' Diel_constant ';
elseif  OptType == 7
    ProName = '   Band_gap    ';
elseif  OptType == 8
    ProName = '   Diel_gap    ';
elseif  OptType == 9
    ProName = '  Mag_moment   ';
elseif  OptType == 10
    ProName = ' Quasientropy  ';
elseif  OptType == 11
    ProName = ' Birefringence ';
elseif  OptType == 14
    ProName = 'Thermoele_const';
elseif  OptType == 1101
    ProName = '  Bulk_Modul   ';
elseif  OptType == 1102
    ProName = '  Shear_Modul  ';
elseif  OptType == 1103
    ProName = '  Young_Modul  ';
elseif  OptType == 1104
    ProName = ' Poisson_ratio ';
elseif  OptType == 1105
    ProName = '  Pugh_ratio   ';
elseif  OptType == 1106
    ProName = 'Vicker_hardness';
elseif  OptType == 1107
    ProName = '   Toughness   ';
elseif  OptType == 1108
    ProName = '   Debye_temp  ';
elseif  OptType == 1109
    ProName = 'Sound_velocity ';
elseif  OptType == 1110
    ProName = 'S_wave_velocity';
elseif  OptType == 1111
    ProName = 'P_wave_velocity';
else
    ProName = '    Fitness    ';
end

