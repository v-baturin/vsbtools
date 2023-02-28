function NEB_plotE()

%--------------------------------------------------------
global ORG_STRUC
global POP_STRUC
%--------     System Parameters      ------%
numImages = ORG_STRUC.numImages;

%------------------------------------------
% Not plot E for only 1 Image
if numImages==1
    return;
end
%------------------------------------------%

E=zeros(numImages,1);
for i = 1:numImages
    E(i)= POP_STRUC.POPULATION(i).Enthalpy;
end
%--------------------------------------------------------

Nmatrix=numImages;
FIG=zeros(Nmatrix+1,Nmatrix*3);
%--------------------------------------------------------
MinE=min(E);
MaxE=max(E);
dE=(MaxE-MinE)/(numImages);

Erescale=round( (E-MinE)/dE );

for i=1:numImages
    LINE = Nmatrix+1-Erescale(i);
    FIG(LINE,(i-1)*3+3)=1;
end
FIG(:,2)=2;
%--------------------------------------------------------

fpath=[POP_STRUC.resFolder,'/BarrierFig'];
fp=fopen(fpath,'w');

%fprintf(fp,'VCNEB Calcualtion Step %4d\n\n', POP_STRUC.step);
if 1==ORG_STRUC.CalcType
    fprintf(fp,'VCNEB Calcualtion Step STEPNUM\n\n');
else
    fprintf(fp,'Relaxation Calcualtion Step STEPNUM\n\n');
end
fprintf(fp,'-Figure-:  Energy Barrier \n ^\n');


for i=1:Nmatrix+1
    fprintf(fp,'%1d',FIG(i,:));
    fprintf(fp,'\n');
end
fclose(fp);
unixCmd( ['sed -i.backup "s/0/ /g" ',POP_STRUC.resFolder,'/BarrierFig'] );
unixCmd( ['sed -i.backup "s/1/+/g" ',POP_STRUC.resFolder,'/BarrierFig'] );
unixCmd( ['sed -i.backup "s/2/|/g" ',POP_STRUC.resFolder,'/BarrierFig'] );

fp=fopen(fpath,'a');
for i=1:numImages
    fprintf(fp,'%3s','--o');
end
fprintf(fp,'-->\n');
for i=1:numImages
    fprintf(fp,'%3d',i);
end
fprintf(fp,'\n\n');

fprintf(fp,'      dE (Max_E-Min_E) = %6.4f eV\n',MaxE-MinE);
fprintf(fp,'Image-%2d = %6.4f eV\n',1,E(1));
fprintf(fp,'Image-%2d = %6.4f eV\n',numImages,E(numImages));

fclose(fp);


unixCmd(['sed -i.backup "s/STEPNUM/' num2str(POP_STRUC.step) '/g" ' POP_STRUC.resFolder,'/BarrierFig']);
delete([POP_STRUC.resFolder '/BarrierFig.backup']);
%---------------------------------------------------------
%--     FUNCTION END
