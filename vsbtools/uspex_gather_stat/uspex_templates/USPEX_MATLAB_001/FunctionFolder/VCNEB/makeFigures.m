function makeFigures()


%---------------------------------------------------%
global ORG_STRUC
global POP_STRUC


if 1~=ORG_STRUC.CalcType
    return;
end

if size(ver('Octave'),1)
    OctaveMode = 1;
else
    OctaveMode = 0;
end

numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);

dimension  = ORG_STRUC.dimension;
numDimension = 3*(sumIons+dimension);

%---------------------------------------------------%
cd(ORG_STRUC.homePath)
unixCmd(['rm -f ' ORG_STRUC.resFolder '/BarrierFigure.tif']);
unixCmd(['rm -f ' ORG_STRUC.resFolder '/SelectedEnergyBarrier.tif']);
unixCmd(['rm -f ' ORG_STRUC.resFolder '/EnergyBarrier_vs_Step.tif']);


x=1:numImages;
y=zeros(1,numImages);
for i = 1:numImages
    y(i)=POP_STRUC.POPULATION(i).Enthalpy;
end
h = figure;
plot(x,y,['--bo'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8);
xlabel('Image number');
ylabel('EnergyBarrier (eV)');
title(['VCNEB Calculaiton Step ' num2str(POP_STRUC.step)]);
print(h,'-dtiff','-r300',[ORG_STRUC.resFolder '/EnergyBarrier.tif']);
print(h,'-dpdf' ,[ORG_STRUC.resFolder '/EnergyBarrier.pdf']);
%----------------------------------------------------%


handle = fopen([ORG_STRUC.resFolder '/AuxiliaryFiles/enthalpy_all.dat']);

plotStep=[0];
plotColor=['kmgrby'];

[nothing, maxStepStr]=unix([' cat ' ORG_STRUC.resFolder '/AuxiliaryFiles/enthalpy_all.dat | wc -l']);
maxStep=str2num(maxStepStr);
energyBarrier=zeros(3,maxStep);

plotStep(1)=maxStep;
for i = [10 50 200]
    plotStep(end+1)=floor((maxStep-i)/10)*10;
end


item = 1;
plotItem=1;
try
    while 1
        tmp = fgetl(handle);
        xy = sscanf(tmp,'%g');
        energyBarrier(1,item)=xy(1);
        energyBarrier(2,item)=max(xy(3:end))-xy(3);
        energyBarrier(3,item)=max(xy(3:end))-xy(end);
        if item==1
            minStep=xy(1);
            plotData(1).step=xy(1);
            plotData(1).x=1:xy(2);
            plotData(1).y=xy(3:end);
            plotItem=2;
        elseif  item==maxStep
            plotData(plotItem).step=xy(1);
            plotData(plotItem).x=1:xy(2);
            plotData(plotItem).y=xy(3:end);
            plotItem=plotItem+1;
        else
            for i = plotStep
                if  xy(1)==i+minStep
                    plotData(plotItem).step=xy(1);
                    plotData(plotItem).x=1:xy(2);
                    plotData(plotItem).y=xy(3:end);
                    plotItem=plotItem+1;
                end
            end
        end
        item = item+1;
        if tmp == -1
            break
            fclose(handle);
        end
    end
catch
end



legendStr=[];
status = fclose(handle);
h = figure;
hold on
for i = 1:length(plotData)
    x=plotData(i).x;
    y=plotData(i).y;
    plot(x,y,['--' plotColor(i) 'o'],'LineWidth',1.5,...
        'MarkerEdgeColor','k',...
        'MarkerFaceColor',[plotColor(i)],...
        'MarkerSize',8);
    legendStr{i}=['step ' num2str(plotData(i).step)];
end
if OctaveMode == 0
    legend(legendStr,'FontSize',12);
else
    legend(legendStr);
end

xlabel('Image Number');
ylabel('EnergyBarrier (eV)');
%print(h,'-dtiff','-r300',[ORG_STRUC.resFolder '/SelectedEnergyBarrier.tif']);
print(h,'-dpdf' ,[ORG_STRUC.resFolder '/SelectedEnergyBarrier.pdf']);

%-------------------------------
h = figure;
hold on

x =energyBarrier(1,[1:floor((end-1)/15)+1:end-1,end]);
y1=energyBarrier(2,[1:floor((end-1)/15)+1:end-1,end]);
y2=energyBarrier(3,[1:floor((end-1)/15)+1:end-1,end]);

plot(x,y1,['--bs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',8);
plot(x,y2,['--gs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8);
if OctaveMode == 0
    legend({'initial->final','final->initial'},'FontSize',15);
else
    legend('initial->final','final->initial');
end

xlabel('VCNEB Step');
ylabel('EnergyBarrier (eV)');
print(h,'-dpdf', [ORG_STRUC.resFolder '/EnergyBarrier_vs_Step.pdf']);


%------------------------------
a=zeros(1,numImages);
b=zeros(1,numImages);
c=zeros(1,numImages);

alpha=zeros(1,numImages);
beta=zeros(1,numImages);
gamma=zeros(1,numImages);
for i = 1:numImages
    Lattice = latConverter(POP_STRUC.POPULATION(i).LATTICE);
    if size(Lattice,1) == 1
        Lattice = Lattice';
    end
    Lattice(4:6) = Lattice(4:6)*180/pi;
    a(i)=Lattice(1);  b(i)=Lattice(2);  c(i)=Lattice(3);
    alpha(i)=Lattice(4);   beta(i)=Lattice(5);   gamma(i)=Lattice(6);
end

h = figure;
subplot(2,1,1); plot(1:numImages,a,['--bs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',8);
hold on
subplot(2,1,1); plot(1:numImages,b,['--rs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',8);
subplot(2,1,1); plot(1:numImages,c,['--gs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8);
title('Image Structral Parameters')
if OctaveMode == 0
    legend({'a','b','c'},'Location','SouthEastOutside');
else
    legend('a','b','c');
end
ylabel('length (Angstrom)');
%

subplot(2,1,2); plot(1:numImages,alpha,['--bs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','b',...
    'MarkerSize',8);
hold on
subplot(2,1,2); plot(1:numImages,beta,['--rs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','r',...
    'MarkerSize',8);
subplot(2,1,2); plot(1:numImages,gamma,['--gs'],'LineWidth',1.5,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor','g',...
    'MarkerSize',8);
if OctaveMode == 0
    legend({'\alpha','\beta','\gamma'},'Location','SouthEastOutside') ;
else
    legend('\alpha','\beta','\gamma') ;
end

xlabel('Image Number');
ylabel('Angle (Degrees)');
print(h,'-dpdf' ,[ORG_STRUC.resFolder '/CellParameters_vs_ImageNumber.pdf']);
