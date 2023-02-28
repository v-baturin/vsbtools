function MINHOP_WriteOUTPUT(Ind_No, resFolder, varcomp, Step, flag)

%To update all the necessary items after each structure is done
%1, OUTPUT.txt
%2, Individuals(_relaxed)
%Qiang Zhu (2014/04/01)

global USPEX_STRUC
  gen    =     USPEX_STRUC.POPULATION(Ind_No).gen;
 symg    =     USPEX_STRUC.POPULATION(Ind_No).symg;
 N_atom  = sum(USPEX_STRUC.POPULATION(Ind_No).numIons);
 enth    =     USPEX_STRUC.POPULATION(Ind_No).Enthalpies(Step)/N_atom;
KPOINTS  =     USPEX_STRUC.POPULATION(Ind_No).K_POINTS(Step,:);
if varcomp == 0
 num     =     USPEX_STRUC.POPULATION(Ind_No).numIons;
else
 num     =     USPEX_STRUC.POPULATION(Ind_No).superCell;
end

if flag == 3 %relaxed
volume   = det(USPEX_STRUC.POPULATION(Ind_No).lat0)/N_atom;
presten  = MattoVec(USPEX_STRUC.POPULATION(Ind_No).PressureTensor0);
else
volume   = det(USPEX_STRUC.POPULATION(Ind_No).LATTICE)/N_atom;
presten  = MattoVec(USPEX_STRUC.POPULATION(Ind_No).PressureTensor);
end

composition = sprintf('%3d',num);
shift=[4, 2, 1]; %so far we only consider 6 component
if size(composition,2)<11
   composition=[composition,blanks(shift(length(num)))];
end



if flag <=2
   fpath =  [resFolder '/OUTPUT.txt'];
   fp = fopen(fpath, 'a+');
   fprintf(fp,     '%4d  [%13s] %14.4f  %14.4f     [%2d %2d %2d]  %3d\n', ...
                  Ind_No, composition, enth, volume, KPOINTS(:), symg);
   fclose(fp);
   if flag == 1
      unixCmd(['echo ' num2str(USPEX_STRUC.POPULATION(Ind_No).Enthalpies,'%10.3f') ' >> ' resFolder '/enthalpies_complete.dat']);
   end
end

if flag == 1
   fpath1 = [resFolder '/Individuals'];
   unixCmd(['echo ' num2str(Ind_No) ' ' num2str(presten,'%10.3f') ' >>' resFolder '/presten.dat']);
elseif flag == 2
   fpath1 = [resFolder '/Individuals_relaxed'];
   unixCmd(['echo ' num2str(Ind_No) ' ' num2str(presten,'%10.3f') ' >>' resFolder '/presten_relaxed.dat']);
elseif flag == 3
   fpath1 = [resFolder '/BESTIndividuals'];
   unixCmd(['echo ' num2str(Ind_No) ' ' num2str(presten,'%10.3f') ' >>' resFolder '/BESTpresten.dat']);
elseif flag == 4
   fpath1 = [resFolder '/BESTIndividuals_relaxed'];
   unixCmd(['echo ' num2str(Ind_No) ' ' num2str(presten,'%10.3f') ' >>' resFolder '/BESTpresten_relaxed.dat']);
end

fp1 = fopen(fpath1, 'a+');
fprintf(fp1,'%3d %4d  [%13s] %14.4f  %14.4f     [%2d %2d %2d]  %3d\n', ...
          gen, Ind_No, composition, enth, volume, KPOINTS(:), symg);
fclose(fp1);
