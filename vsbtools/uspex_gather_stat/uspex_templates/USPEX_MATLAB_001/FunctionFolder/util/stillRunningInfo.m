function stillRunningInfo()


tmp = dir('still_running');
if exist('still_running')
    if (now - tmp.datenum > 1/12)    % 1.00 = 24 hours
        delete('still_running');
    else
        disp('still_runing is present, matlab has to exit, please see NOT_YET file for details');
        unixCmd('echo "This is a warning message"   >> NOT_YET');
        unixCmd('echo "It indicates that you attempted to call matlab"      >> NOT_YET');
        unixCmd('echo "when file still_running is still present"            >> NOT_YET');
        unixCmd('echo "It would be problematic when the numParallel is on"  >> NOT_YET');
        unixCmd('echo "Possible reasons: "                                  >> NOT_YET');
        unixCmd('echo "1, The time interval to call MATLAB is too short"    >> NOT_YET');
        unixCmd('echo "2, MATLAB exits with error"                          >> NOT_YET');
        quit
    end
end

unixCmd('echo "This file is present for two possible reasons"      > still_running');
unixCmd('echo "1, MATLAB is still running"                        >> still_running');
unixCmd('echo "2, MATLAB exits with error"                        >> still_running');
unixCmd('echo "If it stays for long time when numParallel is on"  >> still_running');
unixCmd('echo "Matlab is either in dead loop or exits with error.">> still_running');
unixCmd('echo "please do check it                              "  >> still_running');
