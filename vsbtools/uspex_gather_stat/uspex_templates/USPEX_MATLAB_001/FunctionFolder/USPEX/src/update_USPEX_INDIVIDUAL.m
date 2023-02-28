function update_USPEX_INDIVIDUAL(POPULATION, resFolder, gen, atomType)

%To update all the necessary items after each Individual structure is full relaxed
%Lastly updated by Qiang Zhu (2014/02/18)
global ORG_STRUC
global USPEX_STRUC
global FORCE_STRUC


ID = POPULATION.Number;
USPEX_STRUC.POPULATION(ID).gen     =  gen; %by default we consider it's new struc.
USPEX_STRUC.POPULATION(ID).ToCount =  1;   %by default we consider it's new struc.
USPEX_STRUC.POPULATION(ID).Fitness = 1000; %Initiallization
%Commmon variables
if ORG_STRUC.dimension == 0 && ORG_STRUC.varcomp == 1
    Common_var = {'COORDINATES','LATTICE','numIons','FINGERPRINT','order',...
        'Parents','S_order','symg','K_POINTS','Enthalpies','struc_entr',...
        'howCome','INIT_COORD','INIT_LAT','INIT_FP'};
else
    Common_var = {'COORDINATES','LATTICE','numIons','FINGERPRINT','order',...
        'Parents','S_order','symg','K_POINTS','Enthalpies','struc_entr',...
        'howCome'};
end
for i = 1:length(Common_var)
   eval(['USPEX_STRUC.POPULATION(ID).' Common_var{i} ' = POPULATION.' Common_var{i} ';']);
end

%Privatie variables
Private_var = {'cell', 'numMols', 'MOLECULES', 'numBlocks', 'Surface_numIons',...
           'gap','Fphon','hardness','TE_property', 'mag_moment','magmom_ions',...
           'magmom_ini', 'ldaU','dielectric_tensor', 'birefringence',...
           'elasticMatrix','elasticProperties'};
for i = 1:length(Private_var)
   if isfield(POPULATION, Private_var{i})
      eval(['USPEX_STRUC.POPULATION(ID).' Private_var{i} ' = POPULATION.' Private_var{i} ';']);
   end
end
%some tricky variables
if isfield(POPULATION, 'Vol')
   volume = POPULATION.Vol;
else
   volume = det(POPULATION.LATTICE);
end

USPEX_STRUC.POPULATION(ID).Vol = volume;
USPEX_STRUC.POPULATION(ID).density = calcDensity( POPULATION.numIons, atomType, volume);
if ORG_STRUC.platform>0 % in case of platform = 0 this saving takes too much time
    safesave([resFolder '/USPEX.mat'], USPEX_STRUC);
end

% update the FORCE 
if isfield(POPULATION, 'RelaxStep')
	FORCE_STRUC.POPULATION(ID).atomType     = atomType;
	FORCE_STRUC.POPULATION(ID).numIons      = POPULATION.numIons;
	FORCE_STRUC.POPULATION(ID).RelaxStep    = POPULATION.RelaxStep;

	safesave([resFolder '/FORCE.mat'], FORCE_STRUC);
end
