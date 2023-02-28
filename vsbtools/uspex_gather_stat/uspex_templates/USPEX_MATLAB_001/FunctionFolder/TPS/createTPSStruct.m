function TPS_STRUC  = createTPSStruct(popSize)



%---------------------------------------------------------------------------------------
memStruct_TPS={
    'POPULATION'; 'resFolder'; 'DoneOrder';'genDone'; 'iteration';
    'bodyCount';
    'generation';
    'direction';
    'success';
    'howCome';
    'op';
    'HT';'AHT';'BHT';'lastAHT';'lastBHT';
    'shifter';
    'lastrun';
    'amplitudeA2B'; 'amplitudeB2A'; 'magnitudeA2B'; 'magnitudeB2A';
    'successCounter';
    };

memStruct_POPULATION={
    'Done'; 'ToDo'; 'ToKill'; 'JobID'; 'Folder';'Step';
    'aim'; 'success';           %  MD calc aim phase
    'HT';
    'endOp'; 'minOp'; 'maxOp'; 'endOpStep';'minOpStep'; 'maxOpStep'; 
    'orderParamter';  %  Op data from reading results
    };

TPS_STRUC = [];
for  i = 1:length(memStruct_TPS)
    TPS_STRUC = setfield(TPS_STRUC, memStruct_TPS{i}, {} );
end

POPULATION = [];
for  i = 1:length(memStruct_POPULATION)
    POPULATION = setfield(POPULATION, memStruct_POPULATION{i}, {} );
end
TPS_STRUC(1).POPULATION = POPULATION;


%-----------------------------------------------------------------
TPS_STRUC.resFolder = [];
TPS_STRUC.genDone   = 0;
TPS_STRUC.bodyCount = 0;
TPS_STRUC.iteration = 0;
TPS_STRUC.generation= 0;
TPS_STRUC.DoneOrder = zeros(TPS_STRUC.bodyCount, 1);
TPS_STRUC.direction = '';
TPS_STRUC.op = [];
TPS_STRUC.HT = [];
TPS_STRUC.AHT = [];
TPS_STRUC.BHT = [];
TPS_STRUC.lastAHT = [];
TPS_STRUC.lastBHT = [];
TPS_STRUC.howCome = '';
TPS_STRUC.lastrun.success=0;
TPS_STRUC.lastrun.direction='';

TPS_STRUC.shifter.dH=[];
TPS_STRUC.shifter.T =[];
TPS_STRUC.shifter.randSeed=[];
TPS_STRUC.shifter.accept=[];

TPS_STRUC.successCounter.TPS  = 0;  % all TPS iterations 
TPS_STRUC.successCounter.A2B  = 0;  % every A2B iteration 
TPS_STRUC.successCounter.B2A  = 0;  % every B2A iteration


for i = 1:popSize
    TPS_STRUC.POPULATION(i).Folder= 0;
    TPS_STRUC.POPULATION(i).Done  = 0;
    TPS_STRUC.POPULATION(i).ToDo  = 1;
    TPS_STRUC.POPULATION(i).ToKill= 0;
    TPS_STRUC.POPULATION(i).JobID = 0;  
    TPS_STRUC.POPULATION(i).Step  = 1;
    
    TPS_STRUC.POPULATION(i).minOp = [];
    TPS_STRUC.POPULATION(i).maxOp = [];
    TPS_STRUC.POPULATION(i).endOp = [];
    TPS_STRUC.POPULATION(i).minOpStep = [];
    TPS_STRUC.POPULATION(i).maxOpStep = [];
    TPS_STRUC.POPULATION(i).endOpStep = [];
    
    TPS_STRUC.POPULATION(i).aim   =[];
    TPS_STRUC.POPULATION(i).success=[];
    
    TPS_STRUC.POPULATION(i).orderParamter = [];
end
