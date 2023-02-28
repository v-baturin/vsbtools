function checkTPSSuccess_op()

global ORG_STRUC
global TPS_STRUC


OpACriteria = ORG_STRUC.opCriteria(1);
OpBCriteria = ORG_STRUC.opCriteria(2);

opInd = [0,0];
if strcmp(TPS_STRUC.direction, 'A2B')
    opInd(1) = OpACriteria;
    opInd(2) = OpBCriteria;
else
    opInd(1) = OpBCriteria;
    opInd(2) = OpACriteria;
end

%
%---  orderParameter Format ----
%
%  STEP   Op
%    1   0.1
%

%-- Op reached stauts
for i = 1:length(TPS_STRUC.POPULATION)
    
    if ~isempty( TPS_STRUC.POPULATION(i).orderParamter )
        [minOp,seq] = min( TPS_STRUC.POPULATION(i).orderParamter(:,2) );
        TPS_STRUC.POPULATION(i).minOp       = minOp;
        TPS_STRUC.POPULATION(i).minOpMDStep = TPS_STRUC.POPULATION(i).orderParamter(seq,1);
        
        [maxOp,seq] = max( TPS_STRUC.POPULATION(i).orderParamter(:,2) );
        TPS_STRUC.POPULATION(i).maxOp       = maxOp;
        TPS_STRUC.POPULATION(i).maxOpMDStep = TPS_STRUC.POPULATION(i).orderParamter(seq,1);
        
        TPS_STRUC.POPULATION(i).endOp       = TPS_STRUC.POPULATION(i).orderParamter(end,2);
        TPS_STRUC.POPULATION(i).endOpStep   = TPS_STRUC.POPULATION(i).orderParamter(end,1);
        
        %         if  TPS_STRUC.POPULATION(i).Done
        %             if minOp < minOpCriteria
        %                 TPS_STRUC.POPULATION(i).reachA = 1;
        %             else
        %                 TPS_STRUC.POPULATION(i).reachA = 0;
        %             end
        %             if maxOp < maxOpCriteria
        %                 TPS_STRUC.POPULATION(i).reachB = 1;
        %             else
        %                 TPS_STRUC.POPULATION(i).reachB = 0;
        %             end
        %         end
        
    end
end

%
%-- Succeed or Fail for shoot MD
%

for i = 1:length(TPS_STRUC.POPULATION)
    if TPS_STRUC.POPULATION(i).Done ==1
        %if TPS_STRUC.POPULATION(i).aim == 'A'
        if     abs(opInd(i)- min(opInd)) < 1E-3
            if     TPS_STRUC.POPULATION(i).endOp < opInd(i) %TPS_STRUC.POPULATION(i).minOp < opInd(i)
                TPS_STRUC.POPULATION(i).success = 1;
            elseif  TPS_STRUC.POPULATION(i).maxOp > max(opInd)
                TPS_STRUC.POPULATION(i).success =-1;
            else
                TPS_STRUC.POPULATION(i).success = 0;
            end
        elseif abs(opInd(i)- max(opInd)) < 1E-3
            if      TPS_STRUC.POPULATION(i).endOp > opInd(i) %TPS_STRUC.POPULATION(i).maxOp > opInd(i)
                TPS_STRUC.POPULATION(i).success = 1;
            elseif  TPS_STRUC.POPULATION(i).minOp < min(opInd)
                TPS_STRUC.POPULATION(i).success =-1;
            else
                TPS_STRUC.POPULATION(i).success = 0;
            end
        end
        %end
        
    end
end

%-- Succeed or Fail for the whole MD
if TPS_STRUC.bodyCount == 2
    if      (TPS_STRUC.POPULATION(1).success + TPS_STRUC.POPULATION(2).success) == 2
        TPS_STRUC.success =  1;
    elseif  (TPS_STRUC.POPULATION(1).success + TPS_STRUC.POPULATION(2).success) ==-2
        TPS_STRUC.success = -1;
    else
        TPS_STRUC.success =  0;
    end
    
    TPS_STRUC.op = revertData( TPS_STRUC.POPULATION(1).orderParamter, TPS_STRUC.POPULATION(2).orderParamter );
    TPS_STRUC.HT = revertData( TPS_STRUC.POPULATION(1).HT, TPS_STRUC.POPULATION(2).HT );
end


%-- Update the HT 
for i = 1:length( TPS_STRUC.POPULATION )
    if  TPS_STRUC.POPULATION(i).aim == 'A'
        TPS_STRUC.AHT = TPS_STRUC.POPULATION(i).HT(end,:);
    end
    if  TPS_STRUC.POPULATION(i).aim == 'B'
        TPS_STRUC.BHT = TPS_STRUC.POPULATION(i).HT(end,:);
    end
end




cd (ORG_STRUC.homePath);
safesave ('Current_ORG.mat', ORG_STRUC)
safesave ('Current_TPS.mat', TPS_STRUC)

%===================================%
%
%
%
function allOp = revertData(op1, op2)

allOp = [op1; op2];
len   = size(op1,1);


for i = 1:len
    allOp(i,1)     =-op1(len-i+1,1);
    allOp(i,2:end) = op1(len-i+1,2:end);
end

