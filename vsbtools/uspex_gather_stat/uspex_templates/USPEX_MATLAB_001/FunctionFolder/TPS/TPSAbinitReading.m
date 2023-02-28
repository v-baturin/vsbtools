function Error = TPSAbinitReading(code, ind)


global ORG_STRUC
global TPS_STRUC

Error=0;
%
% Run user defined op calculation commnads
%
if ORG_STRUC.orderParaType   %-- orderParaType=1,user define
    [status, msg]  = unix(ORG_STRUC.cmdCalcOp);
    if status > 0
        disp(msg);
        error('Failed to run user defined op calculation script... USPEX quits.');
        quit();
    else
        try
            TPS_STRUC.POPULATION(ind).orderParamter = importdata(ORG_STRUC.orderParameterFile);
        catch
            error('Error in importing op.dat file.... USPEX quits');
            quit();
        end
    end
else
    if      code==4   % Lammps
        calcfp_lammps=[ORG_STRUC.USPEXPath,'/FunctionFolder/TPS/calcfpDistance_lammps.py'];
        opt1 = ['fpA.dat rdf.dat > fpADistance.dat'];
        opt2 = ['fpB.dat rdf.dat > fpBDistance.dat'];
        [status1, msg]  = python_uspex(calcfp_lammps, opt1);
        [status2, msg]  = python_uspex(calcfp_lammps, opt2);
        
        
        if status1+status2 > 0
            disp(msg);
            error('Failed to run fp calculation script... USPEX quits.');
            quit();
        else
            try
                orderParamterA = importdata('fpADistance.dat');
                orderParamterB = importdata('fpBDistance.dat');
                
                TPS_STRUC.POPULATION(ind).orderParamter = [orderParamterA, orderParamterB(:,2)];
            catch
                error('Error in importing fp.dat file.... USPEX quits');
                quit();
            end
        end
    elseif code==12   % CP2K
                [status, msg]  = unix(ORG_STRUC.cmdCalcfp);
        if status > 0
            disp(msg);
            error('Failed to run user defined op calculation script... USPEX quits.');
            quit();
        else
            try
                TPS_STRUC.POPULATION(ind).orderParamter = importdata(ORG_STRUC.orderParameterFile);
            catch
                error('Error in importing op.dat file.... USPEX quits');
                quit();
            end
        end
    end
end


%
% Import Enthalpy & Tempture data for shifter
%
if ~isempty( ORG_STRUC.cmdCalcHT )
    [status, msg]  = unix(ORG_STRUC.cmdCalcHT);
end
if status > 0
    disp(msg);
    error('Failed to run user defined op calculation script... USPEX quits.');
    quit();
else
    try
        TPS_STRUC.POPULATION(ind).HT = importdata('HT.dat');
    catch
        error('Error in import HT.dat file.... USPEX exit');
        quit();
    end
end
