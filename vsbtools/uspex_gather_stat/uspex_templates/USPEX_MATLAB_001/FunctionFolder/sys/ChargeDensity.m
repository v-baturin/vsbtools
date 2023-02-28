function ChargeDensity(ID, Step)
%%% This file stores charge densities for last step, if ORG_STRUC.chargeDensity = 1
global ORG_STRUC

if ORG_STRUC.chargeDensity == 1
  chgsum=[ORG_STRUC.USPEXPath,'/FunctionFolder/Tool/chgsum.pl'];
  totStep = length(ORG_STRUC.abinitioCode);
  ChgDen = ([ORG_STRUC.homePath,'/' ORG_STRUC.resFolder,'/ChargeDensity']);
  if ~exist(ChgDen)
     mkdir(ChgDen);
  end
  if Step == totStep      %% Only Save Charge densities for last step of relaxation
     unixCmd(['perl ' chgsum ' AECCAR0 AECCAR2']);
     unixCmd(['mv AECCAR0 ' ChgDen, '/CoreChg-' ID '']);
     unixCmd(['mv AECCAR2 ' ChgDen, '/ValChg-' ID '']);
     unixCmd(['mv CHGCAR_sum ' ChgDen, '/TotChg-' ID '']); 
  end
end
