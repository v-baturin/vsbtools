function Finish_PSO()

global POP_STRUC
global ORG_STRUC

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

fprintf(fp,'This calculation runs %4d structural relaxations\n', POP_STRUC.bodyCount);
plotE_PSO(POP_STRUC.resFolder);
fprintf(fp,'\n', datestr(now));
fprintf(fp,'            Job Finished at       %30s\n', datestr(now));
POP_STRUC.generation = ORG_STRUC.numGenerations + 1;
save('Current_POP.mat', 'POP_STRUC')
fclose(fp);
