function checkTPSSuccess()


global ORG_STRUC
global TPS_STRUC

%orderParaType  (Method to calculate the order parameter. 0: fingerPrint, 1. user define  )
if     ORG_STRUC.orderParaType == 1
    checkTPSSuccess_op();
elseif ORG_STRUC.orderParaType == 0
    checkTPSSuccess_fp();
else
    error('We dont know what method to use to run Order Parameter check... USPEX quits.')
end

disp(['<-> Checking TPS MD status ' ]);
disp(['  --> Direction : ' TPS_STRUC.direction ]);


%
% Output the success/failed status
%
if TPS_STRUC.POPULATION(1).success==1
    status='Success';
else
    status='Failed';
end
disp(['   -> Calc-1 Aim: ' TPS_STRUC.direction(1), ', ' status ]);
if TPS_STRUC.POPULATION(2).success==1
    status='Success';
else
    status='Failed';
end 
disp(['   -> Calc-2 Aim: ' TPS_STRUC.direction(3), ', ' status ]) ;

%
% Update the success/failed counter
%

if TPS_STRUC.success==1
    % all TPS iterations
    if TPS_STRUC.successCounter.TPS <= 0
        TPS_STRUC.successCounter.TPS =1;
    else
        TPS_STRUC.successCounter.TPS = TPS_STRUC.successCounter.TPS+1;
    end
    % every A2B/B2A iteration
    if      TPS_STRUC.direction(3) == 'B'
        if TPS_STRUC.successCounter.A2B <= 0
            TPS_STRUC.successCounter.A2B =1;
        else
            TPS_STRUC.successCounter.A2B = TPS_STRUC.successCounter.A2B+1;
        end
    else   %TPS_STRUC.direction(3) == 'A'
        if TPS_STRUC.successCounter.B2A <= 0
            TPS_STRUC.successCounter.B2A =1;
        else
            TPS_STRUC.successCounter.B2A = TPS_STRUC.successCounter.B2A+1;
        end
    end
else % TPS_STRUC.success==0
    % all TPS iterations
    if TPS_STRUC.successCounter.TPS >= 0
        TPS_STRUC.successCounter.TPS =-1;
    else
        TPS_STRUC.successCounter.TPS = TPS_STRUC.successCounter.TPS-1;
    end
    % every A2B/B2A iteration
    if      TPS_STRUC.direction(3) == 'B'
        if TPS_STRUC.successCounter.A2B >= 0
            TPS_STRUC.successCounter.A2B =-1;
        else
            TPS_STRUC.successCounter.A2B = TPS_STRUC.successCounter.A2B-1;
        end
    else   %TPS_STRUC.direction(3) == 'A'
        if TPS_STRUC.successCounter.B2A >= 0
            TPS_STRUC.successCounter.B2A =-1;
        else
            TPS_STRUC.successCounter.B2A = TPS_STRUC.successCounter.B2A-1;
        end
    end 
end

disp(' ')
disp(' --> TPS MD iteration success counter :')
disp(['   -> Continuous TPS Success/Failed : ' num2str(TPS_STRUC.successCounter.TPS )]);
disp(['   -> Continuous A2B Success/Failed : ' num2str(TPS_STRUC.successCounter.A2B )]);
disp(['   -> Continuous B2A Success/Failed : ' num2str(TPS_STRUC.successCounter.B2A )]);



%
%  The first traj checking :
%
if (TPS_STRUC.success==0) && (TPS_STRUC.iteration ==0)
    error('The first trajectory is not correct, pls check your MD restart file or the order parameter data set ...')
end
