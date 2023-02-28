function runShooter()



global ORG_STRUC

disp(' ')
disp('<-> Run Shoot here :')


shooterPy_mole   = [ORG_STRUC.USPEXPath,'/FunctionFolder/TPS/shooter_mole.py'];
shooterPy_atomic = [ORG_STRUC.USPEXPath,'/FunctionFolder/TPS/shooter_atomic.py'];


%-- write the files for python input
tmpFolder = [ORG_STRUC.homePath '/CalcFoldTemp' ];
cd(tmpFolder);


% write InputShoot File
writeInputShoot(); 

%-- copy the MD restart file
unixCmd([' rm ' ORG_STRUC.restartFile ]); % delete the existing MD restart file.
restartFile = [ORG_STRUC.homePath '/' ORG_STRUC.resFolder '/' ORG_STRUC.restartFile];
copyMDRestartFile(restartFile, 'ascii.data');

%
%-- Different shooter test
%
atomic_success = 0;
mol_success    = 0;

if atomic_success ==0
    try
        python_uspex( shooterPy_atomic, [' input_shoot > shoot.log'] );
        atomic_success = 1;
    catch
        atomic_success = 0;
    end
end

if atomic_success ==0
    try
        python_uspex( shooterPy_mole,   [' input_shoot > shoot.log'] );
        mol_success = 1;
    catch
        mol_success = 0;
    end
end

if mol_success + atomic_success ==0
   error('Error to run USPEX TPS shooter, please check your MD restart file... USPEX quits.') 
end

disp('  -> See CalcFoldTemp/shoot.log for more details.')
disp('<-> Shoot DONE.')

