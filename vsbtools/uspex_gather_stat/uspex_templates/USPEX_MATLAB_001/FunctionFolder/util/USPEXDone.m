function writeUSPEXDone(calcType)

        if exist('still_running')
            delete('still_running');
            delete('POSCAR');
            delete('POSCAR_order');
        end
        disp('USPEX IS DONE!');
        unixCmd('echo 1 > USPEX_IS_DONE');
        disp(' ');
        quit
