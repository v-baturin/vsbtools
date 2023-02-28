function checkTPSSuccess_fp()

global ORG_STRUC
global TPS_STRUC


fpACriteria = ORG_STRUC.opCriteria(1);
fpBCriteria = ORG_STRUC.opCriteria(2);


%
%---  orderParameter Format for finger print ----
%
%  STEP   fpA   fpB
%    1    0.99  0.99 
%

%-- fp reached stauts
for i = 1:length(TPS_STRUC.POPULATION)
    
    if ~isempty( TPS_STRUC.POPULATION(i).orderParamter )
        [minfpA,seqA] = min( TPS_STRUC.POPULATION(i).orderParamter(:,2) );
        [minfpB,seqB] = min( TPS_STRUC.POPULATION(i).orderParamter(:,3) );
        TPS_STRUC.POPULATION(i).minOp       = [minfpA,minfpB] ;
        TPS_STRUC.POPULATION(i).minOpMDStep = [TPS_STRUC.POPULATION(i).orderParamter(seqA,1), TPS_STRUC.POPULATION(i).orderParamter(seqB,1),];
        
        [maxfpA,seqA] = max( TPS_STRUC.POPULATION(i).orderParamter(:,2) );
        [maxfpB,seqB] = max( TPS_STRUC.POPULATION(i).orderParamter(:,3) );
        TPS_STRUC.POPULATION(i).maxOp       = [maxfpA,maxfpB];
        TPS_STRUC.POPULATION(i).maxOpMDStep = [TPS_STRUC.POPULATION(i).orderParamter(seqA,1), TPS_STRUC.POPULATION(i).orderParamter(seqB,1),];
        
        TPS_STRUC.POPULATION(i).endOp       = TPS_STRUC.POPULATION(i).orderParamter(end,2:3);
        TPS_STRUC.POPULATION(i).endOpStep   = TPS_STRUC.POPULATION(i).orderParamter(end,1);

    end
end

%
%-- Succeed or Fail for shoot MD
%

for i = 1:length(TPS_STRUC.POPULATION)
    if TPS_STRUC.POPULATION(i).Done ==1
        if TPS_STRUC.POPULATION(i).aim == 'A'
            if     TPS_STRUC.POPULATION(i).endOp(1)  > fpACriteria
                TPS_STRUC.POPULATION(i).success = 1;
            elseif  TPS_STRUC.POPULATION(i).maxOp(1) > fpACriteria
                TPS_STRUC.POPULATION(i).success =-1;
            else
                TPS_STRUC.POPULATION(i).success = 0;
            end
        else 
            if      TPS_STRUC.POPULATION(i).endOp(2) > fpBCriteria
                TPS_STRUC.POPULATION(i).success = 1;
            elseif  TPS_STRUC.POPULATION(i).maxOp(2) > fpBCriteria
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

