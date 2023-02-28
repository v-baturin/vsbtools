function target = Read_Gaussian(flag, ID)
if     flag == 0
        [nothing, results] = unix(' grep "Normal termination of Gaussian" Gauss_output');
        [nothing, results1] = unix('./getStuff  Gauss_output " Maximum Force   " 6 | tail -1');
        [nothing, results2] = unix('./getStuff  Gauss_output " RMS     Force   " 6 | tail -1');
        [nothing, results3] = unix('./getStuff  Gauss_output " Maximum Displacement " 6 | tail -1');
        [nothing, results4] = unix('./getStuff  Gauss_output " RMS     Displacement " 6 | tail -1');
        if ~isempty(results)
		if (length(results1) == 4) && (length(results2) == 4) && (length(results3) == 4) && (length(results4) == 4)
                	disp(['| Gaussian completely Done at ' ID]);
                	target = 1;
	        else
        	        disp('Gaussian is not completely Done');
                	unix(['cp Gauss_output ERROR-OUTPUT-' ID]);
               		target = 0;
		end
        else
                disp('Gaussian is not completely Done');
                unix(['cp Gauss_output ERROR-OUTPUT-' ID]);
                target = 0;
        end
elseif flag == -1
        [nothing, results1] = unix('grep "I DONT THE ERROR" Gauss_output');
        if ~isempty(results1)
                unix(['cp Gauss_output ERROR-OUTPUT-' ID]);
                target = 0;
        else
                disp(['* Gaussian completely Done at ' ID]);
                target = 1;
        end
elseif flag ==  1
        [nothing, results] = unix('./getStuff Gauss_output " SCF Done:" 6 | tail -1');
        target = str2num(results);
%	1 a.u =   27.211396 eV
	target = target*27.211396; %eV
end

