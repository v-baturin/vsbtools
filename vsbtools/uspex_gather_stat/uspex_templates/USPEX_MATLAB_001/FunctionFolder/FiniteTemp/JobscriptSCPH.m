function JobscriptSCPH(NDISP)
% prepare runscph script for SCAILD run  
%TTD

global POP_STRUC;
global ORG_STRUC;

loop = ORG_STRUC.SCPHParam.SCPHloop;

fo = fopen('runscph.sh','w');
fprintf(fo,'#!/bin/sh\n');
fprintf(fo,'#PBS -N scphjob\n');
fprintf(fo,'#PBS -l nodes=1:ppn=24\n'); 
fprintf(fo,'#PBS -l walltime=240:00:00\n'); 
fprintf(fo,'#PBS -q low\n');
fprintf(fo,'#PBS -j oe\n');
fprintf(fo,'cd $PBS_O_WORKDIR\n');
fprintf(fo,'cp INCAR_6 INCAR\n');
%number of displacement 
for ii = 1:NDISP
	posc = ['POSCAR_START' num2str(ii)];
	fprintf(fo,'%s %s %s\n','cp',posc,'POSCAR');
	fprintf(fo,'vasp >out\n');
	fprintf(fo,'./extract_force\n');
	if ii ==1 
		fprintf(fo,'cp FORCESI FORCES.0\n');
		fprintf(fo,'cp FORCESI KRAFTER.1\n');
	else
		krafter = ['KRAFTER.' num2str(ii)];
	       fprintf(fo,'%s %s\n','cp FORCESI',krafter);	
	end 
end 
fprintf(fo,'cp POSCAR_REF POSCAR\n');
fprintf(fo,'scph\n');
fprintf(fo,'cp POSCARTEMP1 POSCARTEMP1.0\n');
fprintf(fo,'cp FREQ FREQ.0\n');
iteration = [1:loop];
fprintf(fo,'%s %s\n','for i in' ,num2str(iteration));
fprintf(fo,'do\n');
fprintf(fo,'cp POSCARTEMP1 POSCAR \n');
fprintf(fo,'rm WAVECAR\n');
fprintf(fo,'vasp >out\n');
fprintf(fo,'./extract_force\n');
fprintf(fo,'cp FORCESI KRAFTER \n');
fprintf(fo,'cp FORCESI FORCES.$i \n');
fprintf(fo,'cp POSCAR_REF POSCAR\n');
fprintf(fo,'scph\n');
fprintf(fo,'cp MEMORYFILE MEMORYFILE.$i\n');
fprintf(fo,'cp POSCARTEMP1 POSCARTEMP1.$i \n');
fprintf(fo,'cp OUTCAR OUTCAR.$i \n');
fprintf(fo,'cp FREQ FREQ.$i \n');
fprintf(fo,'done \n');
fprintf(fo,'cp POSCAR_SCPH CONTCAR\n');
fprintf(fo,'cp OSZICAR_old OSZICAR\n');
fprintf(fo,'rm -rf OUTCAR* WAVECAR\n');
fprintf(fo,'touch SCPH_IS_DONE\n');
fclose(fo);
