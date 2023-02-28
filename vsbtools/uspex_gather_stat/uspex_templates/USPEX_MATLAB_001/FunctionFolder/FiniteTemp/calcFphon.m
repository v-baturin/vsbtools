% This function returns the phonon free energy calculated using SCPH code 
%TTD 

function [Fphon] = calcFphon()
 
 
 
% detemine the current folder name
  currentdir = pwd;
  [upperpath,dir] = fileparts(currentdir);
 

%get Free energy from scph calculation
[temp,results] = unix('./getFphon');

results = str2num(results);
Fphon = results(end);
 
