function     Write_LAMMPS_TPS(Ind_No)



global TPS_STRUC
global ORG_STRUC


aim = TPS_STRUC.POPULATION(Ind_No).aim;

LAMMPS_clean();
unixCmd(['cp -rf ../Specific/* . ']);

%
%
%
maxMDLoop  = ORG_STRUC.maxMDLoop;
randMDLOOP = ceil(maxMDLoop*rand(1));

[status1, msg] = unix(['cat lammps.in_aim' aim ' > lammps.in']);
[status2, msg] = unix(['sed -i.backup "s/TPSMaxLOOP/',  num2str(maxMDLoop+2),  '/g"  lammps.in']);
[status3, msg] = unix(['sed -i.backup "s/TPSRandLOOP/', num2str(randMDLOOP),   '/g"  lammps.in']);

if status1+status2+status3> 0
    error = ['param file is not present for lammps.in_aim' TPS_STRUC.POPULATION(Ind_No).aim];
    disp(msg)
    disp(error)
    save ([ ORG_STRUC.resFolder '/ERROR_param.txt'],'error')
    quit
end


restartFile = [ORG_STRUC.homePath '/CalcFoldTemp/ascii.new'];
copyMDRestartFile(restartFile);


%
%
%
function LAMMPS_clean()


unixCmd('rm lammps*');
unixCmd('rm *dat');
unixCmd('rm *log*');
unixCmd('rm *msd*');
unixCmd('rm *crn*');
unixCmd('rm *backup*');




