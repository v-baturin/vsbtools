function update_USPEX_INDIVIDUAL(POPULATION, resFolder, gen, atomType, flag)

%To update all the necessary items after each Individual structure is full relaxed
%Lastly updated by Qiang Zhu (2014/02/18)
global USPEX_STRUC

ID = POPULATION.Number;
if flag == 1
USPEX_STRUC.POPULATION(ID).gen         =  gen; %by default we consider it is new struc.
USPEX_STRUC.POPULATION(ID).ToCount     =  0;   %by default we consider it is new struc.
elseif flag == 2
USPEX_STRUC.POPULATION(ID).ToCount     =  1;   %by default we consider it is new struc.
end
%Commmon variables
Common_var = {'COORDINATES','LATTICE','numIons','FINGERPRINT','order',...
'PressureTensor', 'Parents', 'symg','K_POINTS','Enthalpies', 'superCell',...
'coords0', 'lat0', 'PressureTensor0'};
for i = 1:length(Common_var)
   eval(['USPEX_STRUC.POPULATION(ID).' Common_var{i} ' = POPULATION.' Common_var{i} ';']);
end

safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
