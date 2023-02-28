function MINHOP_Finish()
%add a complete sign, generally work for every type of calculation
%Last updated by Qiang Zhu (2013/10/17)

global POP_STRUC
global ORG_STRUC

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');

fprintf(fp,'This calculation runs %4d structural relaxations\n', POP_STRUC.bodyCount);
fprintf(fp,'\n', datestr(now));
fprintf(fp,'            Job Finished at       %30s\n', datestr(now));

fclose(fp);

POP_STRUC.generation = ORG_STRUC.numGenerations + 1;
safesave('Current_POP.mat', POP_STRUC)
disp(' ')
delete('still_running');
