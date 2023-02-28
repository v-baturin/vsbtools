function WriteGenerationBackup()

global POP_STRUC
global ORG_STRUC
global USPEX_STRUC
global POOL_STRUC
global CLUSTERS

safesave ([POP_STRUC.current_dir '/POP_STRUC.mat'],  POP_STRUC)
safesave ([POP_STRUC.current_dir '/ORG_STRUC.mat'],  ORG_STRUC)
safesave ([POP_STRUC.current_dir '/USPEX.mat'],    USPEX_STRUC)
safesave ([POP_STRUC.current_dir '/POOL.mat'],      POOL_STRUC)

if (ORG_STRUC.dimension==0) && (ORG_STRUC.varcomp==1) %%for varcomp 001
    safesave ([POP_STRUC.current_dir '/CLUSTERS.mat'], CLUSTERS)
end
