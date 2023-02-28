%check wether the calculation is done and extract results 
%Tekalign
function [target] = Read_SCPH(flag, ID)
%Success 1 ; fail 0
% flag:
%   0 - check whether the calculation is finished
%   3 - get Fphon

    global ORG_STRUC;
    global POP_STRUC;

    currentdir = pwd;
    [upperpath,directory] = fileparts(currentdir);

    if flag == 0

        if exist('SCPH_IS_DONE') && exist('CONVERGENCE')

            target = 1;
        else
            target = 0;
        end

elseif flag == 3 % get the free energy from SCAILD calculation
        target = calcFphon();
end
