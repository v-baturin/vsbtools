function     Read_LAMMPS_TPS(ind)



global TPS_STRUC
global ORG_STRUC


% run user defined op calculation commnads

if ORG_STRUC.orderParaType   %-- orderParaType=1,user define
    [status, msg]  = unix(cmdOrderParameter);
    if status > 0
        disp(msg);
        disp('Failed to run user defined op calculation script... USPEX quits.');
    else
        TPS_STRUC.POPULATION(ind).orderParamter = importdata('op.dat');
    end
else
    % user the finger print to run op.
    % More codes needed ...
    
    TPS_STRUC.POPULATION(ind).orderParamter = cosdistantce();
end





