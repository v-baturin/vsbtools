function WriteTPSOutput()


global ORG_STRUC
global TPS_STRUC



success='FF';
if TPS_STRUC.POPULATION(1).success==1
    success(1) = 'T';
else
    success(1) = 'F';
end
if TPS_STRUC.POPULATION(2).success==1
    success(2) = 'T';
else
    success(2) = 'F';
end


HowCome=  TPS_STRUC.howCome;
iteration = TPS_STRUC.iteration;
aim    =  [TPS_STRUC.POPULATION(1).aim, TPS_STRUC.POPULATION(2).aim];


if TPS_STRUC.shifter.accept == 1
    accept  ='Y';
else
    accept  ='N';
end
if TPS_STRUC.success==1
    dH      = TPS_STRUC.shifter.dH;
    tempT   = TPS_STRUC.shifter.T;
end

ampA2B = TPS_STRUC.amplitudeA2B;
ampB2A = TPS_STRUC.amplitudeB2A;

if ORG_STRUC.orderParaType == 0 % fp
    if TPS_STRUC.POPULATION(1).aim=='A'
        op1 = TPS_STRUC.POPULATION(1).endOp(1);
        op2 = TPS_STRUC.POPULATION(2).endOp(2) ;       
    else
        op1 = TPS_STRUC.POPULATION(1).endOp(2);
        op2 = TPS_STRUC.POPULATION(2).endOp(1);
    end
else
    op1 = TPS_STRUC.POPULATION(1).endOp(1);
    op2 = TPS_STRUC.POPULATION(2).endOp(1);
end
%
%  Write OUTPUT file
%

fpath = [ORG_STRUC.homePath,'/',ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a');


%--OUTPUT format
%Iteration Calc1(Suc)->Calc2(Suc)   Amplitude    Shootor  Shifter:   dH  -  Temp
%                                (A2B)    (B2A)           (accept)(Kcal/mol) (K)

if      TPS_STRUC.success==1 && TPS_STRUC.shifter.accept==1
    fprintf(fp, '>');
elseif  TPS_STRUC.success==1 && TPS_STRUC.shifter.accept==0
    fprintf(fp, '.');
else
    fprintf(fp, ' ');
end
fprintf(fp, ' %5d       %s(%5.4f,%s)->%s(%5.4f,%s)     %7s %7.4g %7.4g',...
    iteration, aim(1), op1, success(1), aim(2), op2, success(2), HowCome, ampA2B, ampB2A );

if TPS_STRUC.success == 1
    fprintf(fp, '      %s  %10.5g %5.2f  \n',accept,dH,tempT);
else
    fprintf(fp, '\n');
end

fclose(fp);


