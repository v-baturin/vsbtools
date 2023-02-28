function Write_DFTB(Ind_No)

global POP_STRUC
global ORG_STRUC

numIons     = POP_STRUC.POPULATION(Ind_No).numIons;
Step        = POP_STRUC.POPULATION(Ind_No).Step;
COORDINATES = POP_STRUC.POPULATION(Ind_No).COORDINATES;
LATTICE     = POP_STRUC.POPULATION(Ind_No).LATTICE;
atomType    = ORG_STRUC.atomType;

Write_DFTB_gen(atomType, numIons, LATTICE, COORDINATES);

try
    unix(['cat dftb_' num2str(Step) ' > dftb_in.hsd']);
catch
    disp(['dftb is not present for step ' num2str(Step)]);
    quit
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KPOINTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fp = fopen('dftb_in.hsd','a+');
fprintf(fp, '  MaxAngularMomentum = {\n');
for i = 1:length(atomType)
    if atomType(i)==1
       fprintf(fp, '   H = "s"\n');
    elseif atomType(i)==5
       fprintf(fp, '   B = "p"\n');
    elseif atomType(i)==6
       fprintf(fp, '   C = "p"\n');
    elseif atomType(i)==7
       fprintf(fp, '   N = "p"\n');
    elseif atomType(i)==8
       fprintf(fp, '   O = "p"\n');
    elseif atomType(i)==9
       fprintf(fp, '   F = "p"\n');
    elseif atomType(i)==15
       fprintf(fp, '   P = "d"\n');
    elseif atomType(i)==16
       fprintf(fp, '   S = "p"\n');
    elseif atomType(i)==17
       fprintf(fp, '  Cl = "p"\n');
    end
end
fprintf(fp, '  }\n');

[Kpoints, Error] = Kgrid(LATTICE, ORG_STRUC.Kresol(Step), ORG_STRUC.dimension);
fprintf(fp, '  KPointsAndWeights = SupercellFolding {\n');
fprintf(fp, '  %2d  0  0\n', Kpoints(1));
fprintf(fp, '   0 %2d  0\n', Kpoints(2));
fprintf(fp, '   0  0 %2d\n', Kpoints(3));
fprintf(fp, '   0  0  0 \n');
fprintf(fp, '  }\n');
fprintf(fp, '}\n');

fclose(fp);
