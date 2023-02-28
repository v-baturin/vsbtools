function mlFilterTrain()

global ORG_STRUC;
global POP_STRUC;
    
assert(ORG_STRUC.mlPeerFilter > 1, ...
       'Valid values for mlPeerFilter are 2, 3.');

assert(ORG_STRUC.mlPeerFilter < 4, ...
       'Valid values for mlPeerFilter are 2, 3.');

% TODO: impl. model 'knr' and 'nn'
model_names = {'knr', 'krr', 'svr', 'nn'};
model = model_names{ORG_STRUC.mlPeerFilter};

res_folder = [ORG_STRUC.homePath '/' POP_STRUC.resFolder];
uspexml_folder = [res_folder '/' 'uspexml'];

model_file = ['uspexml-' model '-gen' num2str(POP_STRUC.generation) '.pkl'];
model_path = [uspexml_folder '/' model_file];
uspex_mat = [res_folder '/' 'USPEX.mat'];
disp(['Training a ' model ' model with structures before generation ' ...
		    num2str(POP_STRUC.generation) '.']);
train_cmd = ['uspexml train -m ' model ' -f ' model_path ...
             ' -d ' uspex_mat ' -v 5 -i 128 -l 8'];
             % ' -p ' ORG_STRUC.resFolder ' -v 5 -i 128 -l 8'];
[s, w] = unix(train_cmd);
if 0 ~= s
  disp(['Error running ' train_cmd]);
end
    
