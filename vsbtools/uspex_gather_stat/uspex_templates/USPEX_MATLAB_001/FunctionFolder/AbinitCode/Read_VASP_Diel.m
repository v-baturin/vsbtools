function Target = Read_VASP_Diel(Type)

% reads dielectric susceptibility tensor from OUTCAR file. Format:
% MACROSCOPIC STATIC DIELECTRIC TENSOR (including local field effects in DFT)
% -------------------------------------
%           4.061     0.000     0.000
%           0.000     4.061     0.000
%           0.000     0.000     4.061
% -------------------------------------
% MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC CONTRIBUTION
% -------------------------------------
%          13.471     0.000     0.000
%           0.000    13.471     0.000
%           0.000     0.000    13.462
% -------------------------------------
global ORG_STRUC

Opt_Sign = ORG_STRUC.opt_sign;

% d_s = 0;
ErrorSB = 0;  % 0 = no Error; 1 = Has a Erorr
ErrorST = 0;  % 0 = no Error; 1 = Has a Erorr
E_Check = 0;  % 0 = no Tensor; 1 = Tensor found
I_Check = 0;  % 0 = no Tensor; 1 = Tensor found
Birefringence = zeros(1,6);
d_s = zeros(1,6);
d_s_ion = zeros(1,6);
[fileHandle, msg] = fopen('OUTCAR');
%%%-------------
try
    while 1
        tline = fgetl(fileHandle);
        if ~ischar(tline)
            break
        end
        if ~isempty(strfind(tline, 'MACROSCOPIC STATIC DIELECTRIC TENSOR ('))
            tline = fgetl(fileHandle);
            tline1 = fgetl(fileHandle);
            tline2 = fgetl(fileHandle);
            tline3 = fgetl(fileHandle);
            e_tensor = str2num(vertcat(tline1, tline2, tline3));
            ElectroTens = e_tensor;
            d_s(1) = e_tensor(1,1);
            d_s(2) = e_tensor(2,2);
            d_s(3) = e_tensor(3,3);
            d_s(4) = e_tensor(1,2);
            d_s(5) = e_tensor(1,3);
            d_s(6) = e_tensor(2,3);
            DiagMat = eig(e_tensor);
            Birefringence(1) = sqrt(max(DiagMat)) - sqrt(min(DiagMat));
            for a = 1 : length(DiagMat)
                if ( DiagMat(a) > 10 ) | ( DiagMat(a) < 1 ) | (d_s(a) < 0)
                    Birefringence = WrongAnswer(Opt_Sign);
                    ErrorSB = 1;
                    break;
                end
            end
            E_Check = 1;
        end
        
        if ~isempty(strfind(tline, 'MACROSCOPIC STATIC DIELECTRIC TENSOR IONIC'))
            tline = fgetl(fileHandle);
            tline1 = fgetl(fileHandle);
            tline2 = fgetl(fileHandle);
            tline3 = fgetl(fileHandle);
            e_tensor = str2num(vertcat(tline1, tline2, tline3));
            IonicTens = e_tensor;
            d_s_ion(1) = e_tensor(1,1);
            d_s_ion(2) = e_tensor(2,2);
            d_s_ion(3) = e_tensor(3,3);
            d_s_ion(4) = e_tensor(1,2);
            d_s_ion(5) = e_tensor(1,3);
            d_s_ion(6) = e_tensor(2,3);
            I_Check = 1 ;
        end
    end
    status = fclose(fileHandle);
catch
    ErrorST = 1;
end

if E_Check == 0
    ErrorSB = 1;
elseif I_Check == 0
    ErrorST = 1;
else
    TotalTens = ElectroTens + IonicTens;
    TotDiagMat = eig(TotalTens);
    for a = 1 : length(TotDiagMat)
        if ( d_s_ion(a) < 0) | (TotDiagMat(a) < 1 ) | (TotDiagMat(a) > 99)
            ErrorST = 1;
        end
    end
    
    try
        [tmp, imaginary_f] = unix('./getStuff OUTCAR f/i 4'); % check for imaginary frequencies
        if ~isempty(imaginary_f)
            i_f = abs(str2num(imaginary_f));
            if (i_f(end) > 0.3) | (i_f(end-1) > 0.3) | (i_f(end-2) > 0.3)
                ErrorST = 1;
            end
        end
    catch
        ErrorST = 1;
    end
end

if Type == 'B'
    if ErrorSB == 0
        Target = Birefringence;
    else
        Target = WrongAnswer(Opt_Sign);
    end
elseif Type == 'D'
    if (ErrorSB + ErrorST) == 0
        Target = d_s + d_s_ion;
    else
        Target = WrongAnswer(Opt_Sign);
    end
end


function d_s = WrongAnswer(optSign)
d_s(1:6) = optSign * -100000 * [ 1 1 1 1 1 1 ];


