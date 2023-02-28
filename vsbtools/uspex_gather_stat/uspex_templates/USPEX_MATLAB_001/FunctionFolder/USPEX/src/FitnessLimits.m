function fitnessRange = FitnessLimits(optType)

fitnessRange = [-99999 99999];

if abs(optType) == 1
    fitnessRange = [-999999999 99999];
elseif abs(optType) == 2
    fitnessRange = [ 0  100];
elseif abs(optType) == 3 
    fitnessRange = [-200 3.4];
elseif abs(optType) == 4
    fitnessRange = [-100  0];
elseif abs(optType) == 5
    if optType < 0
        fitnessRange = [ 0 9999];
    else
        fitnessRange = [-9999 0];
    end
elseif abs(optType) == 6 
    if optType < 0
        fitnessRange = [ 1    20];
    else
        fitnessRange = [-9999 -1];
    end
elseif abs(optType) == 7
    fitnessRange = [ -15  0 ];
elseif abs(optType) == 8
    fitnessRange = [ -999 0 ];
elseif abs(optType) == 9
    fitnessRange = [ -100 0 ];
elseif abs(optType) == 10
    fitnessRange = [ -100 0 ];
elseif abs(optType) == 11
    if optType < 0
        fitnessRange = [ 0  10 ];
    else
        fitnessRange = [ -10 0 ];
    end
elseif abs(optType) == 14
    if optType < 0
        fitnessRange = [ 0  10];
    else
        fitnessRange = [-100 0];
    end
elseif abs(optType) == 1101
    fitnessRange = [ -9999 1 ];
elseif abs(optType) == 1102
    fitnessRange = [ -9999 0 ];
elseif abs(optType) == 1103
    fitnessRange = [ -9999 0 ];
elseif abs(optType) == 1104
    fitnessRange = [ -1  0.5 ];
elseif abs(optType) == 1105
    fitnessRange = [ -10   0 ];
elseif abs(optType) == 1106
    fitnessRange = [ -200  4 ];
elseif abs(optType) == 1107
    fitnessRange = [ -100  0 ];    
elseif abs(optType) == 1108
    fitnessRange = [ -9999 0 ];
elseif abs(optType) == 1109
    fitnessRange = [-99999 0 ];
elseif abs(optType) == 1110
    fitnessRange = [-99999 0 ];
elseif abs(optType) == 1111
    fitnessRange = [-99999 0 ];
end


