function Write_SCPH(Ind_No)
% prepare  inputs files for SCAILD calculation 
%
%Tekalign T. Debela 

global ORG_STRUC;
global POP_STRUC;
Step	 =  POP_STRUC.POPULATION(Ind_No).Step;
numI  =  POP_STRUC.POPULATION(Ind_No).numIons;
symg     =  POP_STRUC.POPULATION(Ind_No).symg;
count    =  POP_STRUC.POPULATION(Ind_No).Number;
atomType =  ORG_STRUC.atomType;
mass     =  elementMass(atomType);
temp     =  ORG_STRUC.SCPHParam.Temp;
displ    = ORG_STRUC.SCPHParam.DISP;
loop     = ORG_STRUC.SCPHParam.SCPHloop;
MPgrid   = ORG_STRUC.SCPHParam.MPgrid;
dosinputs =ORG_STRUC.SCPHParam.Dosinputs;

%determine current directory

currentdir = pwd;
[upperpath,dir] = fileparts(currentdir);

%keep previous files 

unixCmd('cp CONTCAR_old CONTCAR');%starting structure for scph calculation
%unixCmd('cp OSZICAR_old OSZICAR');
unixCmd('chmod +x collect_Fe clean_scph extract_force');
unixCmd('cp CONTCAR_old POSCAR_SCPH');

casename = [dir '-' num2str(Ind_No)];
%FiniteTempRelax()
%update the coords and lattice after relaxation  
[cords, latt] = Read_VASP_Structure();

%get_primitive()

%generate Kpoints 
%get i,j,k for supercel generation 
[Kpoints, Error] = Kgrid(latt, ORG_STRUC.Kresol(Step), ORG_STRUC.dimension);

if Error == 1  % This LATTICE is extremely wrong, let's skip it from now
	POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
else
       	POP_STRUC.POPULATION(Ind_No).K_POINTS(Step,:)=Kpoints;
end

%generate a supercell
[lats,coords,numIons] = createSupercell(latt,cords,numI,Kpoints(1,1),Kpoints(1,2),Kpoints(1,3));
content = POSCARContent(atomType, count, symg, numIons,lats,coords);
%------------------------------------------------
%unixCmd('cp SPOSCAR POSCAR_REF');

%write KPOINTS
fp = fopen('KPOINTS', 'w');
fprintf(fp,'SCPH\n');
fprintf(fp,'0\n');
fprintf(fp,'M\n');
fprintf(fp,'%4d %4d %4d\n', Kpoints(1,:));
fclose(fp);
%------------------------------------------------
% input for phon to displace atoms 
%
fi1 = fopen('INPHON','w');
fprintf(fi1,'LCENTRAL = .F.\n');
fprintf(fi1,'LSUPER = .T.\n');
fprintf(fi1,'%s %4d\n','NTYPES =',length(atomType)); 
fprintf(fi1,'%s %8.3f\n',' MASS =',mass);
fprintf(fi1,'%s %4d %4d %4d\n', 'NDIM =', Kpoints(1,:));
fprintf(fi1,'%s %4d\n',' DISP = ',displ);
fprintf(fi1,'%s %4d\n',' NPOINTS =',MPgrid);
fprintf(fi1,'%s %4d %s %s %4d %s%s %4d\n','QA = ',MPgrid,';','QB = ',MPgrid,';','QC =',MPgrid);
fprintf(fi1,'%s %4d\n','DOSIN = ',dosinputs(1));
fprintf(fi1,'%s %4d\n','DOSEND = ',dosinputs(2));
fprintf(fi1,'%s %8.4f\n','DOSSTEP =',dosinputs(3));
fprintf(fi1,'%s %8.4f\n','DOSSMEAR = ',dosinputs(4));
fprintf(fi1,'LFORCEOUT =.T.\n');
fclose(fi1);
%get displacement and the ID of displaced atom
unixCmd('cp POSCAR_SCPH POSCAR');
unixCmd('phon');
unixCmd('cp SPOSCAR POSCAR_REF');
unixCmd('rm INPHON');%no loger needed 

%write input file  for scph run
fid2 = fopen('INPHON','w');
fprintf(fid2,'%s %4d\n','NTYPES =',length(atomType));
fprintf(fid2,'%s %8.3f\n','MASS =',mass);
fprintf(fid2,'LSUPER = .F.\n');
fprintf(fid2,'%s %4d %4d %4d\n',' NDIM = ',Kpoints(1,:));
fprintf(fid2,'%s %4d\n', 'DISP = ',displ);
fprintf(fid2,'LFREE = .F.\n');
fprintf(fid2,'%s %4d\n',' TEMPERATURE =',temp);
fprintf(fid2,'LRECIP = .T.\n');
fprintf(fid2,'%s %4d\n',' ND =',max(Kpoints(1,:)));
fprintf(fid2,'%s %4d\n','NPOINTS =',MPgrid);
fprintf(fid2,'LGAMMA = .TRUE.\n');
fprintf(fid2,'%s %4d %s %s %4d %s%s %4d\n','QA = ',-1*MPgrid,';','QB = ',MPgrid,';','QC =',MPgrid);
fprintf(fid2,'%s %4d\n','DOSIN = ',dosinputs(1));
fprintf(fid2,'%s %4d\n','DOSEND = ',dosinputs(2));
fprintf(fid2,'%s %8.4f\n','DOSSTEP =',dosinputs(3));
fprintf(fid2,'%s %8.4f\n','DOSSMEAR = ',dosinputs(4));
fprintf(fid2,'LFORCEOUT =.T.\n');
fclose(fid2);
%prepare new poscars with displaced atoms 
Makedisp(atomType,count, symg, numIons, lats, coords);

