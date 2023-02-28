function runShifter()

global ORG_STRUC
global TPS_STRUC

aut = 2.4188843265E-17;   % second
aul = 5.2917720859E-11;   % meter
auT = 3.1577464E+5;       % kelvin
aue = 4.35974417E-18;     % J
kb1 = 1.3806488E-23;      % J.K-1
kb2 = 1.987204118E-3;     % Kcal.mol-1.K-1
avg = 6.0221415E+23;       % Avogadro number
aum = 9.10938291E-31;     % Kg


aim = TPS_STRUC.direction(end);
if aim=='A'
    currentHT = TPS_STRUC.AHT;
    lastHT    = TPS_STRUC.lastAHT;
else
    currentHT = TPS_STRUC.BHT;
    lastHT    = TPS_STRUC.lastBHT;
end

accept = 0;
if TPS_STRUC.success == 1 
    disp(' ');
    disp('<-> The new trajectory succeeded, run Shift here :');
    %
    %HT  : step, enthalpy(Kcal.mol-1.K-1), temperature
    %
    enthalp     = currentHT(2);
    oldEnth     =    lastHT(2);
    temperature = currentHT(3);
    
    kbt = temperature*kb2; % unit convert here!
    dH  = enthalp-oldEnth;
    
    TPS_STRUC.shifter.dH       = dH;
    TPS_STRUC.shifter.T        = temperature;
    TPS_STRUC.shifter.randSeed = rand(1);
    
    
    disp(['  -> dEnthalpy  = ' num2str(dH) ' (Kcal/mol)']);
    disp(['  -> Temperature= ' num2str(temperature) ' (Kelvin)']);
    
    if enthalp < oldEnth
        disp('  -> New trajectory has lower enthalp, accetped');
        accept = 1;
    elseif rand(1) > ORG_STRUC.shiftRatio
        disp(['  -> New trajectory has been accetped (shiftRatio = ' num2str(ORG_STRUC.shiftRatio) ')']);
        accept = 1;
    else
        %   Metropolis MC step
        disp(['  ->    -dH/(kB*T)= ' num2str(-dH/kbt)]);
        disp(['  -> RandomNumber = ' num2str(TPS_STRUC.shifter.randSeed) ]);
        if exp( -dH/kbt ) > TPS_STRUC.shifter.randSeed
            disp('  -> exp( -dH/T) >  RandomNumber ');
            disp('  -> Metropolis MC step: new trajectory has higher enthalp, still accetped');
            accept = 1;
        else
            disp('  -> exp( -dH/T) <  RandomNumber ');
            disp('  -> Metropolis MC step: new trajectory has higher enthalp, rejected');
            accept = 0;
        end
    end
    disp('<-> Shift Done.');
else
    disp(' ');
    disp('<-> The new trajectory failed, no Shift.');
end

TPS_STRUC.shifter.accept = accept;

