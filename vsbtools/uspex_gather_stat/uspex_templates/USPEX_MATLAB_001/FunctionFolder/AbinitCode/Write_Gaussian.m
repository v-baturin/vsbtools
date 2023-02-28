function Write_Gaussian(Ind_No)
global POP_STRUC
global ORG_STRUC
numIons = POP_STRUC.POPULATION(Ind_No).numIons;
step = POP_STRUC.POPULATION(Ind_No).Step;
unix(['cat gaussinput_' num2str(step) ' > input.gjf']);
direct_coord = POP_STRUC.POPULATION(Ind_No).COORDINATES*POP_STRUC.POPULATION(Ind_No).LATTICE;
fp = fopen('input.gjf','a+');
coordLoop = 1;
if ORG_STRUC.dimension == 0
        for i = 1 : length(numIons)
                for j = 1 : numIons(i)
                        fprintf(fp, '%1s   %12.6f %12.6f %12.6f \n' , megaDoof(ORG_STRUC.atomType(i)), direct_coord(coordLoop,:));
                        coordLoop = coordLoop + 1;
                end
        end
else
        unix(['echo    "This code is interfaced only for clusters"']);
end
fclose(fp);
unix(['echo " " >> input.gjf']);
unix(['cat gaussoption_' num2str(step) ' >> input.gjf']);
