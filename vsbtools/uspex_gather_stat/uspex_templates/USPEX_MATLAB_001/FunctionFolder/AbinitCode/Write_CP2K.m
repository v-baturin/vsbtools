function Write_CP2K(Ind_No)


% USPEX Version 8.3.2
% Change: created

global POP_STRUC
global ORG_STRUC

numIons     = POP_STRUC.POPULATION(Ind_No).numIons;
Step        = POP_STRUC.POPULATION(Ind_No).Step;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;
atomType    = ORG_STRUC.atomType;

try
    unixCmd(['cat cp2k_options_' num2str(Step) ' > cp2k.inp ']);
catch
    disp(['cp2k_options is not present for step ' num2str(Step)]);
    quit
end

Lattice = latConverter(LATTICE);
Lattice = Lattice';
Lattice(4:6) = Lattice(4:6)*180/pi;

unixCmd('echo  \&SUBSYS  > subsys.uspex');
unixCmd('echo     \&CELL >> subsys.uspex');
unixCmd(['echo       ABC [angstrom] ' num2str(Lattice(1:3), 11) ' >> subsys.uspex']);
unixCmd(['echo       ALPHA_BETA_GAMMA [deg] ' num2str(Lattice(4:6), 11) ' >> subsys.uspex']);

if     ORG_STRUC.dimension == 3
       unixCmd(['echo       PERIODIC XYZ >> subsys.uspex ']);
elseif ORG_STRUC.dimension == 2
       unixCmd(['echo       PERIODIC XY >> subsys.uspex ']);
elseif ORG_STRUC.dimension == 0
       unixCmd(['echo       PERIODIC NONE >> subsys.uspex ']);
end

unixCmd('echo     \&END >> subsys.uspex');
unixCmd('echo     \&COORD >> subsys.uspex');
unixCmd('echo        SCALED >> subsys.uspex');

coordLoop = 1;
for i = 1 : length(numIons)
 for j = 1 : numIons(i)
  unixCmd(['echo ' megaDoof(atomType(i)) ' ' num2str(COORDINATES(coordLoop,:), 11) ' >> subsys.uspex']);
  coordLoop = coordLoop + 1;
 end
end

%%  print the Pressure 
unixCmd([' echo EXTERNAL_PRESSURE \[kbar\] ' num2str(ORG_STRUC.ExternalPressure*10) ' > pressure.uspex' ]);
%%

unixCmd('echo     \&END >> subsys.uspex');
unixCmd('echo  \&END SUBSYS >> subsys.uspex');
%unixCmd('echo  \&END FORCE_EVAL >> subsys.uspex');
