function [topology_number, candidate, lat] = random_topology(latVolume,numIons)

global ORG_STRUC

str_latVolume = num2str(latVolume);
sGroup = 'NONE';%spaceGroups(nsym);
str_numIons = num2str(numIons);
sz = num2str(size(numIons,2));
str_Ind_No = '0';
%str_coordinatonNumbers = num2str(coordinatonNumbers);
cd([ORG_STRUC.USPEXPath '/FunctionFolder/topology']);
result_top = python_uspex([ORG_STRUC.USPEXPath '/FunctionFolder/USPEX/src/random_topology.py'],str_Ind_No,str_latVolume,sGroup,sz,str_numIons,1);
cd(ORG_STRUC.homePath);
topology_number = result_top(1);
lat = transpose(reshape(result_top(2:10),3,3));
candidate = transpose(reshape(result_top(11:end),3,[]));
