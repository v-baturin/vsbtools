function Random_M400(Ind_No)
% $Rev$
% $Author$
% $Date$

global POP_STRUC
global ORG_STRUC
global OFF_STRUC

angles_num         = size(POP_STRUC.backbone_atoms, 2);
potentialOffspring = Random_Init_M400(angles_num);

% Fill all necessary variables:
OFF_STRUC.POPULATION(Ind_No).ANGLES   = potentialOffspring;
OFF_STRUC.POPULATION(Ind_No).RESIDUES = POP_STRUC.POPULATION(1).RESIDUES;
OFF_STRUC.POPULATION(Ind_No).Parents = [];
OFF_STRUC.POPULATION(Ind_No).numIons = ORG_STRUC.numIons;
OFF_STRUC.POPULATION(Ind_No).howCome = 'Random';

disp(['Structure ' num2str(Ind_No) ' generated randomly']);

end
