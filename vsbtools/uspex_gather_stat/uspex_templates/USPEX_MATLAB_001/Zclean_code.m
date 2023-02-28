%function Zclean_code()
%This is a utility to clean the matlab files (without indents and comments)

%To list all the directories
disp(' Welcome ');
disp(' The following directories will be searched ');
directory =  {'.', 'FunctionFolder', 'FunctionFolder/sys', 'FunctionFolder/symope', 'FunctionFolder/METADYNAMICS', 'FunctionFolder/PSO', ...
             'FunctionFolder/AbinitCode', 'FunctionFolder/VCNEB', 'FunctionFolder/USPEX', 'FunctionFolder/USPEX/src', ...
       'FunctionFolder/USPEX/300', 'FunctionFolder/USPEX/301', 'FunctionFolder/USPEX/311', ...
       'FunctionFolder/USPEX/310', 'FunctionFolder/USPEX/201', 'FunctionFolder/USPEX/200', ...
     'FunctionFolder/USPEX/M200', 'FunctionFolder/USPEX/M300', 'FunctionFolder/USPEX/000', ...
      'FunctionFolder/USPEX/110', 'FunctionFolder/USPEX/PhaseDiagram'};

% Maltab comments starts with % sign, but we don't remove it if it has the following keywords.
keywords = {'num2str(', 'scanf(', 'printf(', 'disp(', 'strrep(', 'echo', 'strfind(', 'format_header', 'format_data', '_format', ' str(1)~=', 'fprintf (', 'findstr(a, '};

disp(' Searching ............');
disp(' ');

Num_keywords = size(keywords, 2);
Num_directory = size(directory,2);
N_count = 0;

for i=1:Num_directory 
    [nothing, nothing] = unix(['rm -rf ' directory{i} '/.svn']);
    Mfile = dir(directory{i});
    N_Mfile = size(Mfile,1);
    for j=1:N_Mfile
      f_name = Mfile(j).name;
      if size(f_name,2) > 2 %to get rid of . and ..
         if(isempty(findstr(f_name, 'Zclean_code')) & isempty(findstr(f_name, 'submitJob')) & isempty(findstr(f_name, 'checkStatusC'))  ...
			           & isempty(findstr(f_name, 'safesave')) & isempty(findstr(f_name, 'python')) )... 
					   & (f_name(end-1:end)=='.m')  %don't apply it on this file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         f_name = [directory{i} '/' f_name];
         fp = fopen(f_name);
         fp1 = fopen([f_name '-back'],'w');
         while 1
            [tmp, count] = fgets(fp);
            if tmp == -1
               break;
            end;
            tmp1 = '';
            first_space = 1;
            for j = 1 : length(tmp)
               if (tmp(j) == '%')
                  keep = 0;
                  for loop = 1:Num_keywords  %how many keywords
                      if ~isempty(findstr(tmp, keywords{loop}))
                         keep = 1;
                         break;
                      end
                  end
                  if keep == 0
                     tmp1 = [tmp1 tmp(end)];
                     break;
                  end
               end
               if ~isspace(tmp(j))
                  first_space = 0;
               end
               if ~first_space
                  tmp1 = [tmp1 tmp(j)];
               end
            end
         
            count = find(isspace(tmp1));
            if (length(tmp1) > 1) & (length(count) < length(tmp1))
               fwrite(fp1, tmp1);
            end
         end
         fclose(fp);
         fclose(fp1);
         [nothing, nothing] = unix(['mv ' f_name '-back '  f_name]);
         N_count = N_count + 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         end
      end
    end
end
disp('DONE');
disp([num2str(N_count) ' Mfiles have been cleaned up']);
disp(' ');

