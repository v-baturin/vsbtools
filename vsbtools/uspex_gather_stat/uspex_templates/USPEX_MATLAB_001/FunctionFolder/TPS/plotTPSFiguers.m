function plotTPSFiguers()


global ORG_STRUC
global TPS_STRUC


% the MD op


cd([ORG_STRUC.homePath, '/', ORG_STRUC.resFolder]);

op = TPS_STRUC.op;
HT = TPS_STRUC.HT;

%
% The whole figure of the MD Simulation
%
h1 = figure;
set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing

if ORG_STRUC.orderParaType == 1 %user defined op
    subplot(3,1,1); plot(op(:,1),op(:,2),'r-');
    ylabel('Order Parameter');
    title(['TPS iteration-' num2str(TPS_STRUC.iteration), ' (' TPS_STRUC.direction ')'] )

else
    subplot(3,1,1); 
    hold on
    plot(op(:,1),op(:,2),'r-');
    plot(op(:,1),op(:,3),'g-');
    ylabel('Structure Similarity');
    legend('Ref: phase A','Ref: phase B', 'Location','SouthEast');
    title(['TPS iteration-' num2str(TPS_STRUC.iteration), ' (' TPS_STRUC.direction ')'] )

    box on
end

subplot(3,1,2); plot(HT(:,1),HT(:,2),'b-');
ylabel('Free Energy (Kcal/mole)');

subplot(3,1,3); plot(HT(:,1),HT(:,3),'c-');
xlabel('MD steps');
ylabel('Temperature (K)');
print(h1,'-dpdf' , 'thermoTPS.pdf');
print(h1,'-dpdf' , ['Iterations/iter-' num2str(TPS_STRUC.iteration),'/thermoTPS.pdf']);



%
% Detial figure of the MD Simulation near the restart point
%
h2 = figure;
set(gcf,'Visible','off');   % Use this switcher to prevent Matlab foregroung printing

restartStep = find( op(:,1)==0 );
dStep       = round( length(op(:,1))/10 );
step1 = restartStep -dStep;
step2 = restartStep +dStep-1;

if ORG_STRUC.orderParaType == 1 %user defined op
    subplot(3,1,1); plot(op(step1:step2,1),op(step1:step2,2),'r-');
    ylabel('Order Parameter');
    title(['TPS iteration-' num2str(TPS_STRUC.iteration), ' (' TPS_STRUC.direction ')'] )

else
    subplot(3,1,1); 
    hold on
    plot(op(step1:step2,1),op(step1:step2,2),'r-');
    plot(op(step1:step2,1),op(step1:step2,3),'g-');
    ylabel('Structure Similarity');
    legend('Ref: phase A','Ref: phase B', 'Location','SouthEast');
    title(['TPS iteration-' num2str(TPS_STRUC.iteration), ' (' TPS_STRUC.direction ')'] )

    box on
end

restartStep = find( HT(:,1)==1 );
dStep       = round( length(HT(:,1))/10 );
step1 = restartStep -dStep;
step2 = restartStep +dStep-1;

subplot(3,1,2); plot(HT(step1:step2,1),HT(step1:step2,2),'b-');
ylabel('Free Energy (Kcal/mole)');


subplot(3,1,3); plot(HT(step1:step2,1),HT(step1:step2,3),'c-');
xlabel('MD steps');
ylabel('Temperature (K)');
print(h2,'-dpdf' , 'thermoTPSDetail.pdf');
print(h2,'-dpdf' , ['Iterations/iter-' num2str(TPS_STRUC.iteration),'/thermoTPSDetail.pdf']);
