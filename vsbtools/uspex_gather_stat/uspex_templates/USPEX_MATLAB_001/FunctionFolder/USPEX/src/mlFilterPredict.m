function [POP] = mlFilterPredict(POPULATION, numFiltered)
  % Parameters:
  %   POPULATION: the structures of the new generation
  %   numFiltered: number of structures left after filtering
  % Return values:
  %   POP: the filtered structures

global ORG_STRUC;
global POP_STRUC;

assert(ORG_STRUC.mlPeerFilter > 1, ...
       'Valid values for mlPeerFilter are 2, 3.');

assert(ORG_STRUC.mlPeerFilter < 4, ...
       'Valid values for mlPeerFilter are 2, 3.');

uspexml_folder = [ORG_STRUC.homePath '/' POP_STRUC.resFolder '/' 'uspexml'];

% TODO: impl. model 'knr' and 'nn'
model_names = {'knr', 'krr', 'svr', 'nn'};
model = model_names{ORG_STRUC.mlPeerFilter};

model_file = ['uspexml-' model '-gen' num2str(POP_STRUC.generation) '.pkl'];
model_path = [uspexml_folder '/' model_file];

POP_PEER = struct();
for whichInd = 1:length(POPULATION)
  if isempty(POPULATION(whichInd).FINGERPRINT)
    LATTICE = POPULATION(whichInd).LATTICE;
    numIons = sum(POPULATION(whichInd).numIons); %!Grey Fp for varcomp
    COORDINATES = POPULATION(whichInd).COORDINATES;
    [Ni, V, dist_matrix, typ_i, typ_j] = makeMatrices(LATTICE, COORDINATES, numIons, 1); %Grey Fp for varcomp
    [order, FINGERPRINT, atom_fing] = fingerprint_calc(Ni, V, dist_matrix, typ_i, typ_j, numIons); %Grey Fp For varcomp
    POP_PEER.POPULATION(whichInd).order =  order;
    POP_PEER.POPULATION(whichInd).FINGERPRINT = FINGERPRINT;
  end
end

disp(['Filtering generated structures using trained ' model ' model.']);

peer_mat = ['POP_PEER' num2str(POP_STRUC.generation) '.mat'];
pop_peer_mat = [uspexml_folder '/' peer_mat];
safesave(pop_peer_mat, POP_PEER)
    
predict_cmd = ['uspexml predict -f ' model_path ' -d ' pop_peer_mat];
    
[ret, text] = unix(predict_cmd);
if 0 == ret
  order = str2num(text);
  count = 1;
  for ind = 1:length(POPULATION)
    if order(ind) <= numFiltered && count <= numFiltered
      POP(count) = POPULATION(ind);
      count = count + 1;
    end
  end
else
  POP = POP_PEER;
  disp(['Error in executing ' predict_cmd]);
end
