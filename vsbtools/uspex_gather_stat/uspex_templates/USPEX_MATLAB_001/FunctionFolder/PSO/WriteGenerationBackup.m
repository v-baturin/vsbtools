function WriteGenerationBackup()

global POP_STRUC
global ORG_STRUC
global PSO_STRUC

safesave ([POP_STRUC.current_dir '/POP_STRUC.mat'],  POP_STRUC)
safesave ([POP_STRUC.current_dir '/ORG_STRUC.mat'],  ORG_STRUC)
safesave ([POP_STRUC.current_dir '/PSO.mat'],        PSO_STRUC)
