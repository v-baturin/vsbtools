function target = Read_DFTB(flag, ID)

target = 0;
if     flag == -1
   [nothing, results] = unix('grep --text "Geometry converged" detailed.out');
   if isempty(results)
      disp('DFTB is not completely Done');
      unixCmd(['cp DFTB_output ERROR-OUTPUT-' ID]);
   else
      target = 1;
   end
elseif flag == 0
   [nothing, results] = unix('./getStuff detailed.out "Total energy" 6');
   [nothing, N_line] = unix('cat geo_end.gen |wc -l');
   [nothing, N_atom] = unix('grep --text F geo_end.gen');
   if (str2num(N_line) > 1) && isempty(findstr(N_atom, 'No such file'))
      N_atom=str2num(N_atom(1:7));
      if N_atom > str2num(N_line)
         disp('geo_en.gen is broken');
      else
         if isempty(results) 
            disp('detail.out in DFTB is not done correctly!!!'); 
            unixCmd(['cp DFTB_output ERROR-OUTPUT-' ID]);
         elseif ~exist('DFTB_output')
            disp('DFTB Output does not exist!');
         else
            target = 1;
         end
      end   
   else
      disp('geo_en.gen is broken');
   end
elseif flag ==  1 
   [nothing, results6] = unix('./getStuff detailed.out "Total energy" 6');
   [nothing, results5] = unix('./getStuff detailed.out "Total energy" 5');
   target = str2num(results6);
   if isempty(target)
      target = str2num(results5);
   end
end

