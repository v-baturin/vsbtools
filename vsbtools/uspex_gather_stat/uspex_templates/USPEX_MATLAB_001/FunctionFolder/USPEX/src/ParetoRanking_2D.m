function [Stru_Num, Distance, ParetoFront] = ParetoRanking_2D(f1, f2, St_Num)
% Main_F
%     structure    f1       f2    pareto-front   Dist  Ranking
%      number                                    
%   39.0000    0.0801         0    1.0000    0.0801    1.0000
%   28.0000    0.0322    0.1747    1.0000    0.1776    2.0000
%   44.0000         0    0.4406    1.0000    0.4406    3.0000
%   24.0000    0.0577    0.3402    2.0000    0.3451    4.0000
%   26.0000    0.0436    0.4573    2.0000    0.4594    5.0000
%   36.0000    0.0278    0.5243    2.0000    0.5250    6.0000
%   38.0000    0.0073    0.5434    2.0000    0.5434    7.0000
%   43.0000    0.0072    0.6447    2.0000    0.6447    8.0000
%   16.0000         0    0.7131    2.0000    0.7131    9.0000
%   27.0000    0.0578    0.3634    3.0000    0.3679   10.0000
%   37.0000    0.0379    0.5390    3.0000    0.5404   11.0000
%   30.0000    0.0185    0.6983    3.0000    0.6985   12.0000
%    1.0000         0    0.7467    3.0000    0.7467   13.0000
%    2.0000    0.1874    0.3756    4.0000    0.4198   14.0000
%   40.0000    0.1101    0.5326    4.0000    0.5439   15.0000
%   25.0000    0.0928    0.5442    4.0000    0.5520   16.0000
%   41.0000    0.0511    0.5797    4.0000    0.5819   17.0000
%   17.0000    0.0366    0.7643    4.0000    0.7652   18.0000
%    8.0000    0.0023    0.7838    4.0000    0.7838   19.0000
%    7.0000         0    0.7998    4.0000    0.7998   20.0000
%   20.0000    0.3044    0.4687    5.0000    0.5589   21.0000
%   12.0000    0.1334    0.6503    5.0000    0.6639   22.0000
%   42.0000    0.1147    0.6646    5.0000    0.6744   23.0000
%   10.0000    0.5138    0.4636    5.0000    0.6920   24.0000
%   32.0000    0.0749    0.7237    5.0000    0.7276   25.0000
%   13.0000    0.0373    0.7660    5.0000    0.7669   26.0000
%   14.0000    0.0276    0.7977    5.0000    0.7982   27.0000
%    4.0000    0.3955    0.4721    6.0000    0.6158   28.0000
%   15.0000    0.1485    0.6570    6.0000    0.6735   29.0000
%    9.0000    0.1366    0.6648    6.0000    0.6787   30.0000
%    3.0000    0.0788    0.7343    6.0000    0.7385   31.0000
%   33.0000    0.0777    0.7895    6.0000    0.7933   32.0000
%   18.0000    0.4192    0.5571    7.0000    0.6972   33.0000
%   19.0000    0.2133    0.6833    7.0000    0.7158   34.0000
%   23.0000    0.3863    0.6763    7.0000    0.7789   35.0000
%   45.0000    0.1733    0.8328    7.0000    0.8507   36.0000
%   29.0000    0.2265    0.9928    8.0000    1.0183   37.0000
%   31.0000    0.2141    1.0000    8.0000    1.0227   38.0000
%    5.0000    0.8400    0.7438    8.0000    1.1220   39.0000
%    6.0000    0.9945    0.7542    9.0000    1.2482   40.0000
%   11.0000    1.0000    0.8002   10.0000    1.2808   41.0000
%   35.0000    0.5138    0.4636  100.0000  100.0000   42.0000
%   34.0000    0.0322    0.1747  100.0000  100.0000   43.0000
%   22.0000    0.5138    0.4636  100.0000  100.0000   44.0000
%   21.0000    0.1874    0.3756  100.0000  100.0000   45.0000
%---------------------------------------------------------------

N = length(f1);
% Find the Utopia point as a reference
utopia_f1 = min(f1);
utopia_f2 = min(f2);

% set the utopia point to zero
f1_left = f1 - utopia_f1;
f2_left = f2 - utopia_f2;

% the hell point
hell_f1 = max(f1_left);
hell_f2 = max(f2_left);

%rescale the points within the box from utopia to the hell
% size is determined by the hell point
%f1_left = f1_left ./ hell_f1;
%f2_left = f2_left ./ hell_f2;

main_tmp = zeros(N,6);
main_tmp(:,1) = St_Num;
main_tmp(:,2) = f1_left;
main_tmp(:,3) = f2_left;
%[main_tmp , duplicate] = removeDuplicateSys(main_tmp);
ParetoCounter = 0;
Counter = 0;
Main_F = [];

%% In this section we determine pareto front, first pareto, second pareto and so on
while size(main_tmp,1) > 0
    ParetoCounter = ParetoCounter + 1;
    Begining = Counter+1;
    for a = size(main_tmp,1) : -1 : 1
        f1 = main_tmp(a,2);
	f2 = main_tmp(a,3);
	if firstParetoFront(f1, f2, main_tmp, a) == 1
	    Counter = Counter +1;
	    main_tmp(a,4) = ParetoCounter;
	    new_tmp(Counter,:) = main_tmp(a,:);
	    remove_row(Counter) = a;
	end
    end
    for a = Begining:Counter
        main_tmp(remove_row(a),:) = [];
    end
end

%%% Measuring Distance of each point from utopia point, 
for a = 1 : size(new_tmp,1)
    new_tmp(a,5) = sqrt((new_tmp(a,2) ^ 2) + (new_tmp(a,3) ^ 2));    % distance
end

%%% Ordering structures according to their pareto front, 
for ParFront = 1 : new_tmp(end,4)    % for a = 1 : last pareto front
    n = 1;
    draft = [];
    for i = 1 : size(new_tmp,1)
	if new_tmp(i,4) == ParFront
	    draft(n,:) = new_tmp(i,:);
	    n = n + 1;
	end
    end
    draf_len = size(draft,1);	
    for i = 1 :draf_len 
	for j = 1 : size(draft,1)-1
	    if draft(j,5) > draft(j + 1,5)
		XC = draft(j,:);
		draft(j,:) = draft(j + 1,:);
		draft(j + 1,:) = XC;
	    end
	end
    end
    Main_len = size(Main_F,1);
    for a = 1 : draf_len
	Main_F(Main_len+a,:) = draft(a,:);
    end 
end
%%% Adding duplicate systems at the end of ranking so they will not be on top( they have already been)
%Main_len = size(Main_F,1);
%dup_len  = size(duplicate,1);
%duplicate(:,4) = 1000;
%duplicate(:,5) = 100;
%for i = 1 : dup_len
%    Main_F(Main_len+i,:) = duplicate(i,:);
%end
Main_F(:,6) = 1:N;

Stru_Num   = Main_F(:,1);
Distance   = Main_F(:,5);
ParetoFront = Main_F(:,4);
%%% Since the some structures are closer to utopia point but they are not in the first pareto front, we
%%% increase their dinstance from utopia point manually to put them in the appropriate position they need to be. 
LEN = size(Distance,1);
for a = 1 : LEN - 1
    if Distance( a + 1 ) < Distance( a )
        DiF = abs(Distance( a ) - Distance( a + 1 ));
        for B = a + 1 : LEN
            Distance(B) = Distance(B) + DiF + 0.003;
        end
    end
end

%%%%%------------------------------------------------
%function [main_tmp , duplicate] = removeDuplicateSys(main_tmp)
%count = 0;
%duplicate = [];
%for a = size(main_tmp,1) : -1 : 2
%    for b = a-1 : -1 : 1
%	if main_tmp(a,2) == main_tmp(b,2) && main_tmp(a,3) == main_tmp(b,3)
%	    count = count + 1;
%	    duplicate(count,:) = main_tmp(a,:);
%	    main_tmp(a,:) = [];
%            break;
%	end
%    end
%end

%%%%%------------------------------------------------
function result = firstParetoFront(f1, f2, main_tmp, ID)
result = 1;
main_tmp(ID,:) = [];
for a = 1 : size(main_tmp,1)
    if (f1 == main_tmp(a,2)) && (f2 == main_tmp(a,3))

    elseif (f1 >= main_tmp(a,2)) && (f2 >= main_tmp(a,3))
        result = 0;
        break;
    end
end

