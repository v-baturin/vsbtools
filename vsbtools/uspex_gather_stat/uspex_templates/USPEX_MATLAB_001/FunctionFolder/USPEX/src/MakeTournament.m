function MakeTournament(howManyProliferate)
% THIS IS IMPORTANT and subject to change. 
% Tournament: the probability table of each struc to be selected as parent 
% the following is quadratic

global ORG_STRUC
global POP_STRUC

tournament = zeros(howManyProliferate,1);
tournament(end) = 1;

for loop = 2 : howManyProliferate
    tournament(end - loop + 1) = tournament(end - loop + 2) + loop ^ 2;
end
ORG_STRUC.tournament = tournament;

% This part was added for pareto, to give all the structures
% of same Pareto front, the same posibility to get chosen!
if ORG_STRUC.paretoRanking ~= 0
    TotalStructures = sum(POP_STRUC.paretoFront);
    Tournament      = zeros(TotalStructures , 1);
    for a = 1 : length(ORG_STRUC.tournament)
	Tournament(a) = ORG_STRUC.tournament(a);
    end
    EndPoint = 0;
    for b = 1 : length(POP_STRUC.paretoFront)
	StartPiont = 2 + EndPoint;  
	EndPoint   = sum(POP_STRUC.paretoFront(1:b));
	Tournament(StartPiont : EndPoint) = 0;	
    end	     
    for a = length(Tournament) : -1 : 1
	if Tournament(a) == 0
	    Tournament(a) = [];
	end
    end
    ORG_STRUC.tournament = Tournament;
end
