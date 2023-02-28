function Makedisp(atomType,count, symg, numIons, lattice, coor)
%prepare POSCAR_STARTs  
%Tekalign T Debela

fid = fopen('DISP','r');%read displaced atomID and displacements 
sizedata = [4 Inf];
displ = fscanf(fid, '%d %f %f %f',sizedata);
displ = displ';
fclose(fid);

NDISP = size(displ,1);
% add the displacement 
for i= 1 : NDISP
	for j = 1 : numIons
       		 if (displ(i,1)==j)
		coor(j,1) = coor(j,1) + displ(i,2);
		coor(j,2) = coor(j,2) + displ(i,3);
		coor(j,3) = coor(j,3) + displ(i,4);

		content = POSCARContent(atomType, count, symg, numIons, lattice,coor);
		filenamep =['POSCAR_START' num2str(i)];
		writeContent2File(filenamep, content,'w');
		end
	end 
end 
%prepare INPUTFLJ input file for SCAILD 
ff = fopen('INPUTFLJ','w');
fprintf(ff,'POT ABI\n');
fprintf(ff,'%s %4d\n','NDISP',NDISP);
fprintf(ff,'DISPL 216\n');

for dd = 1:NDISP
	fprintf(ff,'%4d %8.4f %8.4f %8.4f \n',displ(dd,:));
end
for kk = (3+NDISP):218
	fprintf(ff,'2 0.003 0.0 0.0\n');%%%usless,but should write
end 
fprintf(ff,'RM 3.2323\n');
fprintf(ff,'BETA -1.1363\n');
fprintf(ff,'X1 1.30000\n');
fprintf(ff,'X2 1.65000\n');
fprintf(ff,'B1 0.01331200\n');
fprintf(ff,'B2 -0.24573000\n');
fprintf(ff,'B3 1.9047000\n');
fprintf(ff,'B4 -8.054100\n');
fprintf(ff,'C6 20.21600\n');
fprintf(ff,'C8 -30.8060\n');
fprintf(ff,'C10 28.506000000\n');
fprintf(ff,'EPS 11.285169080\n');
fprintf(ff,'RCUT 30.0\n');
fprintf(ff,'RK  1.00000000  1.05000000  1.55000000  1.60000000  1.65000000 1.70000000  1.7500000\n');
fprintf(ff,'AK -38.4084159  36.9240324 -6.06706900  8.45632500 -4.68932080 7.95571790 -5.6449389\n');
fprintf(ff,'AZ -0.44417136  1.07926027 -0.67936231 -14.902\n');
fprintf(ff,'SAMPLING NGAUSS\n');
fprintf(ff,'NGAUSSAMPLE 1\n');
fprintf(ff,'MAXITTER 400\n');
fprintf(ff,'DSITTER  400\n');
fprintf(ff,'SUPERSAFE .FALSE.\n');
fprintf(ff,'MAXAMP 1.00\n');
fprintf(ff,'SYMETRIZATION .TRUE.\n');
fclose(ff);
JobscriptSCPH(NDISP);
