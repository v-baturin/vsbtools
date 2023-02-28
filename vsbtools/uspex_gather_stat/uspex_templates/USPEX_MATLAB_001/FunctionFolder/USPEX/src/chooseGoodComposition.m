function ind = chooseGoodComposition(tournament, paretoRanking, POP)


count = 0;
goodComposition = 0;

maxTime = length(POP.POPULATION)*10;

while ~goodComposition

    if paretoRanking ~= 0
       ParetoF = find (tournament > randi([0,max(tournament)-1]));
       %%All structures in same paretofront have the same possibility to get chosen!!
       toMutate = randi(POP.paretoFront(ParetoF(end))) + sum(POP.paretoFront(1 : ParetoF(end) - 1));
    else
       toMutate = find (tournament>RandInt(1,1,[0,max(tournament)-1]));
    end
    
    ind = toMutate(end);
    
    numBlocks = POP.POPULATION(ind).numBlocks;
    if ~isempty(numBlocks)  %% for fix composition condition
        goodComposition = CompositionCheck(numBlocks);
    else
        goodComposition = 1;
    end
    
    count = count +1;
    if count > maxTime  
        ind = -1;
        break;
    end
end
