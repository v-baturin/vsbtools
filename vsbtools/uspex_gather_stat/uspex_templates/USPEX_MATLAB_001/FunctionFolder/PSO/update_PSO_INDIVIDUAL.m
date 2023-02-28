function update_PSO_INDIVIDUAL(POPULATION, resFolder, gen, atomType)

%To update all the necessary items after each Individual structure is full relaxed
%Lastly updated by Qiang Zhu (2014/02/18)
global PSO_STRUC

ID = POPULATION.Number;
PSO_STRUC.POPULATION(ID).gen         =  gen; %by default we consider it is new struc.
PSO_STRUC.POPULATION(ID).ToCount     =  1; %by default we consider it is new struc.
PSO_STRUC.POPULATION(ID).Fitness     = 1000; %Initiallization
%Commmon variables
Common_var = {'COORDINATES','LATTICE','numIons','FINGERPRINT','order',...
    'Parents','S_order','symg','K_POINTS','Enthalpies','struc_entr','howCome'};
for i = 1:length(Common_var)
    eval(['PSO_STRUC.POPULATION(ID).' Common_var{i} ' = POPULATION.' Common_var{i} ';']);
end

%Privatie variables
Private_var = {'gap','hardness','mag_moment','dielectric_tensor','elasticMatrix','elasticProperties'};
for i = 1:length(Private_var)
    if isfield(POPULATION, Private_var{i})
        eval(['PSO_STRUC.POPULATION(ID).' Private_var{i} ' = POPULATION.' Private_var{i} ';']);
    end
end
%some tricky variables
if isfield(POPULATION, 'Vol')
    volume = POPULATION.Vol;
else
    volume = det(POPULATION.LATTICE);
end

PSO_STRUC.POPULATION(ID).Vol = volume;
PSO_STRUC.POPULATION(ID).density = calcDensity( POPULATION.numIons, atomType, volume);
safesave([resFolder '/PSO.mat'], PSO_STRUC);
