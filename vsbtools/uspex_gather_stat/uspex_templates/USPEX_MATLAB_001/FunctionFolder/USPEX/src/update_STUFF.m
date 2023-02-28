function OUTPUT = update_STUFF(inputFile)

global ORG_STRUC
global POP_STRUC
global USPEX_STRUC

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Some parameters can be updated from INPUT.txt %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
createORG_Fingerprint(inputFile);
createORG_EA(inputFile);
createORG_Symmetry(inputFile);
createORG_Onthefly(inputFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fraction of each operator can be updated as below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_pro = round((length(POP_STRUC.ranking)) * ORG_STRUC.bestFrac);
if N_pro == 0
    N_pro = 1;
end
howManyProliferate = min([N_pro, length(POP_STRUC.ranking) - POP_STRUC.bad_rank]);

Parent =  {'Random',        'RandTop',  'Heredity', 'H1eredity', ...
           'Permutate', 'AddAtom', 'RemoveAtom',...
           'TransMutate', 'LatMutate',...
           'softmutate',     'Rotate',      'Spin',...
           'SecSwitch', 'ShiftBorder'};
fractions = {'fracRand', 'fracRandTop',  'fracGene', 'fracG1ene',...
             'fracPerm', 'fracAddAtom', 'fracRemAtom',...
             'fracTrans',  'fracLatMut',...
             'fracAtomsMut','fracRotMut','fracSpin',...
             'fracSecSwitch', 'fracShiftBorder'};

minfrac = [0.10, 0.10, 0.05, 0.05, 0, 0.05, 0.05, 0.05, 0, 0.10, 0.10, 0, 0.10, 0.10];
numOperation = length(Parent);
f = zeros(1, numOperation); 

if ~ORG_STRUC.AutoFrac
    for i = 1 : numOperation
        eval([ 'f(1,i) = ORG_STRUC.' fractions{i} ';' ]);
    end
else
    N = zeros(3, numOperation); 
    gen = POP_STRUC.generation;

    if gen > 1
        prev1_gen = USPEX_STRUC.GENERATION(gen-1).ID;  %
        if gen == 2
            prev2_gen = 0;
        else
            prev2_gen = USPEX_STRUC.GENERATION(gen-2).ID;%structures from previous
        end
        for i = 1 : howManyProliferate %for all the good Structures
            C_ID = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).Number;
            accept = 1;
            if C_ID > 0
                for j = prev2_gen + 1 : prev1_gen
                    if ORG_STRUC.dimension == 3 %remove same structs from previous
                        ans = SameStructure_order(j, C_ID, USPEX_STRUC);  
                    else %for low D systems, order is not good
                        ans = SameStructure(j, C_ID, USPEX_STRUC);
                    end
                    if ans
                        accept = 0;
                        break;
                    end
                end
            end
            if accept
                tmp = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).howCome;
                for count = 1 : numOperation
                    if ~isempty(findstr(tmp, Parent{count}))
                        N(1, count) = N(1, count) + 1;
                        break;
                    end
                end
            end
        end

        for i = 1:length(POP_STRUC.POPULATION) %for all the good Structures
            tmp = POP_STRUC.POPULATION(POP_STRUC.ranking(i)).howCome;
            for count = 1 : numOperation
                if ~isempty(findstr(tmp, Parent{count}))
                    N(2, count) = N(2, count) + 1;
                    break;
                end
            end
        end

        count = 0;
        sum1 = 0;
        for i = 1:numOperation
            if N(2, i) > 0
                N(1, i) = max(1, N(1, i)); %make sure N(1,i) >0
                X(i) = N(1, i) / N(2, i);
                count = count + 1;
                sum1 = sum1 + X(i);
            end
        end

        X_mean = sum1/count;
        for i = 1 : numOperation
            if N(2, i) >0
                N(3, i) = round(N(1,i) * (X(i) / X_mean + 1) / 2);
            end
        end
        f = N(3, :) / sum(N(3, :));
    end
    
    for i = 1 : numOperation
        ANS = eval([ 'ORG_STRUC.' fractions{i} '>0' ]);
        eval([ 'f(1,i) = 0.55 * ORG_STRUC.' fractions{i} ' + 0.45 * f(1,i);' ]);
        if ANS > 0
           f(1,i) = max(f(1,i), minfrac(i));
        end
    end
end

f = f / sum(f); %Normalization
for i = 1 : numOperation
    eval( ['ORG_STRUC.' fractions{i} '= f(1,i);'] );
end

for i = 1: numOperation
    table{2*i-1} = Parent{i};
    table{2*i}   = num2str(f(1,i),'%4.2f');
end
OUTPUT = sprintf('    fraction of generation produced by %-14s : %4s\n', table{:});
MakeTournament(howManyProliferate);
