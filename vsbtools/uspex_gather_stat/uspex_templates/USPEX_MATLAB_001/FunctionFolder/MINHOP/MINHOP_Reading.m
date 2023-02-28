function Reading(code, Ind_No, indic)
%This rountine is used to read some necessary tags from the output of Abintio Calculation;
%1, Check is the Calculation correctly done
%2, Read crystal structure (coordinate + lattice)
%3, Read target properties (energy/pressure tensor/dielectric constant, etc)
%4, clean output files in case they will be read at the following cycle
%Update: Qiang Zhu (2013/10/03)

global POP_STRUC
global ORG_STRUC
cd (['CalcFold' num2str(indic)])

TotalStep = length([ORG_STRUC.abinitioCode]);
Step = POP_STRUC.POPULATION(Ind_No).Step;
Gen=POP_STRUC.generation;
ID = ['Gen' num2str(Gen) '-Ind' num2str(Ind_No) '-Step' num2str(Step)]; %For ERROR output
Const_Lat = 1; %fixed lattice, we don't optimize the cell

GoodBad = Read_AbinitCode(code, 0, ID, Step);  %Step 1; Check if the calculation is correctly done
if GoodBad==1
   [POP_STRUC.POPULATION(Ind_No).COORDINATES,POP_STRUC.POPULATION(Ind_No).LATTICE] = Read_Structure(code, Const_Lat);
   POP_STRUC.POPULATION(Ind_No).Enthalpies(Step) = Read_AbinitCode(code, 1, ID, Step);
   POP_STRUC.POPULATION(Ind_No).PressureTensor   = Read_AbinitCode(code, 2, ID, Step);
   POP_STRUC.POPULATION(Ind_No).Step = Step + 1;
   Clean_AbinitCode(code);
else
   if (code~=1) && (code~=8) % If error appears in the empirical code, we just skip this structutre.
      POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 4;
   else
      POP_STRUC.POPULATION(Ind_No).Error = POP_STRUC.POPULATION(Ind_No).Error + 1; %FOR VASP and QE, let's try again
   end
   disp(['PROBLEM_reading Structure' num2str(Ind_No)]);
   disp(pwd);
end
cd ..
