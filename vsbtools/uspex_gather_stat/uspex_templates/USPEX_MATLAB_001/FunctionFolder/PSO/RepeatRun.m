function RepeatRun()

% repeating the same run for statistics

%flag = (POP_STRUC.generation > ORG_STRUC.numGenerations - 0.5);   % flag = 1, otherwise halting criteria won't work and we will get an eternal cycle
flag = 1;
[nothing, rfs] = unix (['./getStuff multiple_runs runs_left 1']);
rfs_1 = str2num(rfs);

if rfs_1 == 1
    quit
end

if (rfs_1 > 1) && (flag)
    rfs_1 = rfs_1 - 1;
    rfs = [num2str(rfs_1) ' : runs_left'];
    unixCmd(['echo "' rfs '" > multiple_runs']);
    unixCmd('rm -r CalcFold*');
    unixCmd('rm -r *.mat');
end

