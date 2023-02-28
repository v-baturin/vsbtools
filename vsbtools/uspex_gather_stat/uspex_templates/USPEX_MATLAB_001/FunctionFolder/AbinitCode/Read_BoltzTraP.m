%This function check if the property calculation succeeds or not 
% Success: 1; No: 0;
% Dong Dong 

% function [target] = checkProperty()
function [target] = Read_BoltzTraP(flag, ID)
% flag:
%   0 - check whether the calculation is finished
%   1 - enthalpy, this should never enter this function
%   2 - power factor, deprecated
%   3 - maximum of tr(ZT)

    global ORG_STRUC;
    global POP_STRUC;

    % detemine the current folder name
    currentdir = pwd;
    [upperpath,directory] = fileparts(currentdir);
  
    % Added by Fei Qi, 2015/07/27
    if flag == 0
        trace_files = dir([directory '-*.trace']);
        if length(trace_files) >= 1
            target = 1;
        else
            target = 0;
        end
    elseif flag == 2
        target = calcThermoelectricProperties();
    elseif flag == 3
        [info, count] = sscanf(ID, 'Gen%d-Ind%d-Step%d');  % decode ID to get Ind_No
        ind_str = num2str(info(2)); % get string of structure index number
        folder = [ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/' 'TEproperties'];
        py_calc_te = [ORG_STRUC.USPEXPath '/FunctionFolder/Tool/calc_te_prop.py'];
        output = python_uspex(py_calc_te, ...
                              ['--output-directory ' folder], ...
                              ['--output-prefix ' ID], ...
                              ['--optimization-goal ' ORG_STRUC.TEparam.goal], ...
                              ['--temperature ' num2str(ORG_STRUC.TEparam.Tinterest)], ...
                              ['--threshold ' num2str(ORG_STRUC.TEparam.threshold)], ...
                              ['--casename ' directory '-' ind_str], currentdir);

        %%% Sometimes output is like " e-05, 0.0481096453231, -0.01749, e-07, 0.0481096453231 "
        %%% In this case str2num(output) = [], and we face trouble, Below we change this kind of
        %%% output to " 1e-05, 0.0481096453231, -0.01749, 1e-07, 0.0481096453231 " to make life
        %%% easier!  ADDeD by Zahed ZX
        counter = 0;
        output_tmp = output;
        for a = 1 : length(output)
            if output(a) == 'e'
                output_tmp(a + counter) = '1';
                output_tmp(end + 1) = ' ';
                output_tmp(a + 1 + counter : end) = output( a : end);
                counter = counter + 1;
            end
        end
        output = output_tmp;
        properties = str2num(output);
        target = properties(end);
    end

  % Commented out Dong Dong's implementation. %% Fei Qi, 2015/07/27
  % if exist([dir '.trace'],'file')
  %    target = 1;
  % else
  %    target = 0;
  % end

