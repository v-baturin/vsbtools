function writeORG()

global ORG_STRUC

fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a+');


% Create a header automatically:
description    = 'Minima Hoping Code for Structure Prediction';
funcFold       = [ORG_STRUC.USPEXPath '/FunctionFolder'];
formatted_rows = createHeader('USPEX_pic.txt', description, funcFold);

for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end


% Cite:
text = {'Please cite the following suggested papers',       ...
    'when you publish the results obtained from USPEX:' ...
    };
formatted_rows = createHeader_wrap(text);
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end

text = {'', ...
    '',  ...
    '' ...
    };
formatted_rows = createHeader_wrap(text, 'left');
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end

text = {datestr(now)};
formatted_rows = createHeader_wrap(text);
for i=1:size(formatted_rows,2)
    fprintf(fp,[formatted_rows{i} '\n']);
end
fprintf(fp,'\n');


fprintf(fp,'            Job Starts at       %30s\n', datestr(now));
fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('Block for system description') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp,'                        Dimensionality  :  %2d\n', ORG_STRUC.dimension);
fprintf(fp,'                        Molecular       :  %2d (1:Yes, 0,No)\n', ORG_STRUC.molecule);
fprintf(fp,'                   Variable Composition :  %2d (1:Yes, 0,No)\n', ORG_STRUC.varcomp);

fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('Block for atomic description') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp,'    There are %1d types of atoms in the system: ', length(ORG_STRUC.atomType));
for i=1:length(ORG_STRUC.atomType)
    fprintf(fp,'%4s', megaDoof(ORG_STRUC.atomType(i)));
end
fprintf(fp,'\n');

fprintf(fp,'    The investigated system is: ');
for i=1:length(ORG_STRUC.atomType)
    fprintf(fp,'%2s_%2d  ', megaDoof(ORG_STRUC.atomType(i)), ORG_STRUC.numIons(i));
end
fprintf(fp,'\n');

fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('Block for evolutionary algorithm') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp,'                 Number of Generations  :    %4d\n', ORG_STRUC.numGenerations);
fprintf(fp,'               General Population Size  :    %4d\n', ORG_STRUC.populationSize);


fprintf(fp, [alignLine('-', 0) '\n']);
fprintf(fp, [alignLine('AB INITIO CALCULATIONS') '\n']);
fprintf(fp, [alignLine('-', 0) '\n']);
total_step = length(ORG_STRUC.abinitioCode);
fprintf(fp,'*  There are %2d local relaxation steps for each individual structures *\n', total_step);
if ORG_STRUC.platform==0
    fprintf(fp,'Step  Abinitio Code                  Execute Command            K-resolution \n');
else
    fprintf(fp,'Step  Abinitio Code    K-resolution \n');
end
for step=1:total_step
    i = ORG_STRUC.abinitioCode(step);
    if i==0
        code=' no relaxation';
    elseif i==1
        code='     VASP     ';
    elseif i==2
        code='    SIESTA    ';
    elseif i==3
        code='     GULP     ';
    elseif i==4
        code='    LAMMPS    ';
    elseif i==5
        code=' Neu Network  ';
    elseif i==6
        code='   DMACRYS    ';
    elseif i==7
        code='    CP2K      ';
    elseif i==8
        code='    PWSCF     ';
    elseif i==9
        code='     ASE      ';
    elseif i==10
        code='     ATK      ';
    elseif i==11
        code='   CASTEP     ';
    else
        fprintf(fp,'***************************************************************\n');
        fprintf(fp,'**      ERROR:    The code you selected is not vaild!        **\n');
        fprintf(fp,'**                 Program STOPS!!!!!!!!!!!                  **\n');
        fprintf(fp,'***************************************************************\n');
    end
    if sum(i==[3 4 5 6])==1 %for GULP, LAMMPS, NN, and DMACRYS
        ORG_STRUC.Kresol(step)=0;
    elseif isempty(ORG_STRUC.Kresol(step))
        fprintf(fp,'********************************************************************\n');
        fprintf(fp,'**      ERROR:   You did not specify Kresolution for each Step    **\n');
        fprintf(fp,'**                 Program STOPS!!!!!!!!!!!                       **\n');
        fprintf(fp,'********************************************************************\n');
    end
    if ORG_STRUC.platform==0
        fprintf(fp,' %2d %12s  %36s  %12.3f\n', step, code, ORG_STRUC.commandExecutable{step}, ORG_STRUC.Kresol(step));
    else
        fprintf(fp,' %2d %12s  %12.3f\n', step, code, ORG_STRUC.Kresol(step));
    end
end
fprintf(fp,'\n');

if ORG_STRUC.platform==0
    fprintf(fp,'The calculations are performed in nonParallel mode on the local machine\n');
elseif ORG_STRUC.platform==1
    fprintf(fp,'The script for job submission is prepared seperately in Submission/*_local.m\n');
elseif ORG_STRUC.platform==2
    fprintf(fp,'The script for job submission is prepared seperately in Submission/*_remote.m\n');
elseif ORG_STRUC.platform==3
    fprintf(fp,'The calculations are performed in Parallel mode on CFN supercomputer\n');
elseif ORG_STRUC.platform==4
    fprintf(fp,'The calculations are performed in Parallel mode on QSH supercomputer\n');
elseif ORG_STRUC.platform==5
    fprintf(fp,'The calculations are performed in Parallel mode on xservDE supercomputer\n');
end
fprintf(fp,'%4d parallel calculations are performed simutaneously\n', ORG_STRUC.numParallelCalcs);
fprintf(fp,'\n');

fclose(fp);
