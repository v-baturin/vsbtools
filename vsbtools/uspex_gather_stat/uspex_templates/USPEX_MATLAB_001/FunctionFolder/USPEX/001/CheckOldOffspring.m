function fp_different = CheckOldOffspring( LATTICE, COORDINATES, numIons, Ind_No, operation )
%CheckOldParents( COORDINATES, LATTICE, numIons, operation ) checks if ...
%brand new structure coincides with one of the parents produced in the past
global USPEX_STRUC
global ORG_STRUC
global OFF_STRUC
global CLUSTERS %#ok<NUSED>

FINGERPRINT = fp_calc_001( LATTICE, COORDINATES, numIons);
toleranceFP = ORG_STRUC.toleranceFing;
fp_different = 1;

if ~strcmp(operation,'AddAtom_001') && ~strcmp(operation,'RemoveAtom_001')
    off_struc_size =  Ind_No-1;
    for past = 1:off_struc_size
        off_numIons = OFF_STRUC.POPULATION(past).numIons;
        off_init_fp = OFF_STRUC.POPULATION(past).INIT_FP;
        
        if isequal(numIons, off_numIons)
            weight = CalcWeight_001(numIons);
            fp_init_dist = cosineDistance(FINGERPRINT, off_init_fp, weight);
            if fp_init_dist < toleranceFP
                eval(['CLUSTERS.repetitions.' operation ' = CLUSTERS.repetitions.'...
                    operation '+1;']);
                %display(['offspring similar to OFF_STRUC structure No ' num2str(past) ' detected']);
                fp_different = 0;
                break;
            end
        end
    end
    
    if fp_different == 1
        uspex_struc_size = length(USPEX_STRUC.POPULATION);
        
        for past = 1:uspex_struc_size
            old_numIons = USPEX_STRUC.POPULATION(past).numIons;
            old_init_fp = USPEX_STRUC.POPULATION(past).INIT_FP;
            old_fp = USPEX_STRUC.POPULATION(past).FINGERPRINT;
            
            if isequal(numIons, old_numIons)
                weight = CalcWeight_001(numIons);
                fp_init_dist = cosineDistance(FINGERPRINT, old_init_fp, weight);
                fp_relax_dist = cosineDistance(FINGERPRINT, old_fp, weight);
                if fp_init_dist < toleranceFP
                    eval(['CLUSTERS.repetitions.' operation ' = CLUSTERS.repetitions.'...
                        operation '+1;']);
                    %display(['offspring similar to INITIAL structure No ' num2str(past) ' detected']);
                    fp_different = 0;
                    break;
                elseif fp_relax_dist < toleranceFP
                    eval(['CLUSTERS.repetitions.' operation ' = CLUSTERS.repetitions.'...
                        operation '+1;']);
                    %display(['offspring similar to RELAXED structure No ' num2str(past) ' detected']);
                    fp_different = 0;
                    break;
                end
            end
        end
    end
end

if fp_different
    OFF_STRUC.POPULATION(Ind_No).INIT_FP = FINGERPRINT;
end
