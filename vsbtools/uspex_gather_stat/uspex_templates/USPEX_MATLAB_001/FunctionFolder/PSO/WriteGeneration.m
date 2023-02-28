function WriteGeneration(resFolder)
%To update all the necessary items at the end of each gen to USPEX_STRUC
%Output the Generation info
%quasiEntropy
%composEntropy
%Lastly updated by Qiang Zhu (2014/02/18)
global PSO_STRUC
fpath = [ resFolder '/OUTPUT.txt'];
fp = fopen(fpath, 'a+');
q_entropy = PSO_STRUC.GENERATION(end).quasiEntropy;

fprintf(fp, [alignLine( sprintf('Quasi entropy = %.4f', q_entropy) ) '\n'] );

if ~isempty(PSO_STRUC.GENERATION(end).composEntropy)
    c_entropy = PSO_STRUC.GENERATION(end).composEntropy;
    fprintf(fp, [alignLine( sprintf('Composition entropy = %.4f', c_entropy) ) '\n'] );
end

fclose(fp);
