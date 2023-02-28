function InitializeVariantsAR()

global ORG_STRUC
global POP_STRUC
global AR_VARIANTS
global CLUSTERS
load CLUSTERS.mat

AR_VARIANTS = [];

covRad = [];
for i = 1 : length(ORG_STRUC.atomType)
    covRad(i) = str2num(covalentRadius(ORG_STRUC.atomType(i)));
end

for i = 1 : length(ORG_STRUC.tournament)
    coef = (length(ORG_STRUC.tournament)-i+1)^2;
    
    same = 0;
    i_pop = POP_STRUC.ranking(i);
    for i_cl = 1 : length(CLUSTERS.Structure)
        if same_structure_001(i_pop, i_cl) %current structure is already in CLUSTERS.Structure database
            same = 1;
            AR_VARIANTS.Structure(i) = CLUSTERS.Structure(i_cl);
            AR_VARIANTS.number(i) = i_cl;

            % recalculation probabilities of atoms according to their ranking numbers
            for j = 1 : sum(AR_VARIANTS.Structure(i).numIons)
                AR_VARIANTS.Structure(i).Atom(j).remove = AR_VARIANTS.Structure(i).Atom(j).remove * coef;
                for k = 1 : length(AR_VARIANTS.Structure(i).Atom(j).Edge)
                    for l = 1 : length(ORG_STRUC.atomType)
                        AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(l) = AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(l) * coef;
                    end
                end
            end
            break;
        end
    end
    
    if same==0 %current structure did not meet before in CLUSTERS.Structure database
        coord       = POP_STRUC.POPULATION(i_pop).COORDINATES;
        numIons     = POP_STRUC.POPULATION(i_pop).numIons;
        lattice     = POP_STRUC.POPULATION(i_pop).LATTICE;
        enthalpy    = POP_STRUC.POPULATION(i_pop).Enthalpies(end);
        fingerprint = POP_STRUC.POPULATION(i_pop).FINGERPRINT;
        
        AR_VARIANTS.Structure(i).coord       = coord;
        AR_VARIANTS.Structure(i).numIons     = numIons;
        AR_VARIANTS.Structure(i).lattice     = lattice;
        AR_VARIANTS.Structure(i).enthalpy    = enthalpy;
        AR_VARIANTS.Structure(i).fingerprint = fingerprint;
        
        CN = coord_numbers(coord, numIons, lattice);
        
        % finding difference between max and current coordination number
        maxCN = [];
        deltaCN = [];
        for ii = 1 : length(numIons)
            begin_atom = sum(numIons(1:ii-1))+1;
            end_atom = sum(numIons(1:ii));
            maxCN(ii) = max(CN(begin_atom : end_atom));
            for jj = begin_atom : end_atom
                deltaCN(jj) = (maxCN(ii) - CN(jj))^2;
            end
        end
        
        for j = 1 : sum(numIons)
            %remove atom
            if sum(deltaCN)>0.0001
                AR_VARIANTS.Structure(i).Atom(j).remove = deltaCN(j)/sum(deltaCN);
            else
                AR_VARIANTS.Structure(i).Atom(j).remove = 1/sum(numIons);
            end
            
            %add atom
            %finding the closest edges to the current atom
            n = 0; %number of the edge
            typeAtom_j = find_atomType(numIons,j);
            coef_r = 1.4;
            while n == 0 %we repeat this procedure until not find at least one atom near the considered atom
                for k = 1 : sum(numIons)
                    dist_cur = distAtoms(coord(k,:)*lattice, coord(j,:)*lattice);
                    typeAtom_k = find_atomType(numIons,k);
                    if (k~=j) && (dist_cur <= coef_r*(covRad(typeAtom_k) + covRad(typeAtom_j)))
                        n = n + 1;
                        AR_VARIANTS.Structure(i).Atom(j).Edge(n).Atom = k;
                    end
                end
                coef_r = coef_r * 1.1;
            end
            
            if sum(deltaCN)>0.0001
                prob = (deltaCN(j)/sum(deltaCN))/(n*length(ORG_STRUC.atomType));
            else
                prob = (1/sum(numIons))/(n*length(ORG_STRUC.atomType));
            end
            for k = 1 : n
                for l = 1 : length(ORG_STRUC.atomType)
                    AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(l) = prob;
                end
            end
        end
        
        % add current structure to CLUSTERS.Structure database
        i_cl = length(CLUSTERS.Structure) + 1;
        CLUSTERS.Structure(i_cl) = AR_VARIANTS.Structure(i);
        AR_VARIANTS.number(i) = i_cl; %link to the CLUSTERS.Structure

        % recalculation probabilities of atoms according to their ranking numbers
        for j = 1 : sum(numIons)
            AR_VARIANTS.Structure(i).Atom(j).remove = AR_VARIANTS.Structure(i).Atom(j).remove * coef;
            for k = 1 : length(AR_VARIANTS.Structure(i).Atom(j).Edge)
                for l = 1 : length(ORG_STRUC.atomType)
                    AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(l) = AR_VARIANTS.Structure(i).Atom(j).Edge(k).add(l) * coef;
                end
            end
        end
    end
end