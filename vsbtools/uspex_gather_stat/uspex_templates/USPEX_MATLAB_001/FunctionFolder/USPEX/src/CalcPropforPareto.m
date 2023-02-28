function Prop =  CalcPropforPareto(pareto)

%%% In this function we gather 3th property 
%%% for Pareto ranking.

global POP_STRUC
global ORG_STRUC

atomType = ORG_STRUC.atomType;

if pareto < 0
opt_sign = -1;
else
opt_sign =  1;
end
pareto = abs(pareto);
Prop = zeros(1, length(POP_STRUC.POPULATION));
for i = 1 : length(POP_STRUC.POPULATION)
    Volume = det(POP_STRUC.POPULATION(i).LATTICE);
    numIons = POP_STRUC.POPULATION(i).numIons;
    if pareto == 2 % volume
        Prop(i) = Volume/sum(numIons);
    elseif pareto == 3 % hardness
        if ~isempty(POP_STRUC.POPULATION(i).hardness)
            Prop(i) = -1*POP_STRUC.POPULATION(i).hardness;
        else
            Prop(i) = 100000;
        end
    elseif pareto == 4 % Structural Order
        if ~isempty(POP_STRUC.POPULATION(i).S_order)
            Prop(i) = -1*POP_STRUC.POPULATION(i).S_order;
        else
            Prop(i) = 100000;
        end
    elseif pareto == 5  % density
        Prop(i) = -1*calcDensity(numIons, atomType, Volume);
    elseif pareto == 6 % dielectric tensor (or rather Tr()/3)
        Prop(i) = -1*sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3;
    elseif pareto == 7 % bandgap
        Prop(i) = -1*POP_STRUC.POPULATION(i).gap;
    elseif pareto == 8 % Tr(dielectric tensor)/3 multiplied by gap^2 value
        Egc = 4; %critical bandgap for semiconductor and insulator, eV
        if POP_STRUC.POPULATION(i).gap >= Egc
            Prop(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^2; %Insulator
        else
            Prop(i) = -1*(sum(POP_STRUC.POPULATION(i).dielectric_tensor(1:3))/3)*(POP_STRUC.POPULATION(i).gap/Egc)^6; %Semiconductor
        end
    elseif pareto == 9
        Prop(i) = -1*POP_STRUC.POPULATION(i).mag_moment/Volume;
    elseif pareto == 10 % structure entropy
        if ~isempty(POP_STRUC.POPULATION(i).struc_entr)
            Prop(i) = -1*POP_STRUC.POPULATION(i).struc_entr;
        else
            Prop(i) = 100000;
        end
    elseif pareto == 11
	Prop(i) = -1*POP_STRUC.POPULATION(i).birefringence; 
    elseif pareto == 14 % thermoelectric property, replaced power factor
        Prop(i) = -1*POP_STRUC.POPULATION(i).TE_property;
    elseif (pareto > 1100) && (pareto < 1112)
        whichPara= mod(pareto,110);
        for i = 1 : length(POP_STRUC.POPULATION)
            if isempty(POP_STRUC.POPULATION(i).elasticProperties) || (POP_STRUC.POPULATION(i).elasticProperties(end)==0)
                Prop(i)=NaN;
            else
                Prop(i) = -1*POP_STRUC.POPULATION(i).elasticProperties(whichPara);
            end
        end
    end            
end

Prop = opt_sign*Prop; % change the mode of optimization (minimization <=> maximization)


