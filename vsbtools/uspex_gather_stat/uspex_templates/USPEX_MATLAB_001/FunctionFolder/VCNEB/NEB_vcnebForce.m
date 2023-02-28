function NEB_vcnebForce()

%--------------------------------------%
global ORG_STRUC
global POP_STRUC

%--------------------------------------%
numImages =     ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);

dimension  = ORG_STRUC.dimension;
numDimension = 3*(sumIons+dimension);

%-----------------------------------------------------------------%
F_pro = zeros(numDimension,1);
F_ela = zeros(numDimension,1);
F_neb = zeros(numDimension,1);
%-----------------------------------------------------------------%

% disp('in NEB_vcnebForce.m')
% ORG_STRUC.pickupImages
%  ORG_STRUC.optMethodCIDI

isCIDIImage = NEB_isTransitionState();
% TS(CI) Images > 0
% LM(DI) Images < 0

ORG_STRUC.whichCI=[];
ORG_STRUC.whichDI=[];

for N=1:numImages
    F = POP_STRUC.POPULATION(N).totCaForceVector;
    tuo = POP_STRUC.POPULATION(N).tuo;
    %disp('Tuo:')
    %disp(tuo)
    %disp(F)
    F_pro = ( F - ( tuo*(tuo'*F) ) ) .*ORG_STRUC.whetherConstraint ;
    %----------------------------------------------------------%
    if (1==N) || (numImages==N)
        k1 = 0;
        k2 = 0;
        F_ela = -tuo*0;
    else
        k1 = 0.5*( POP_STRUC.POPULATION(N).springK + POP_STRUC.POPULATION(N+1).springK );
        k2 = 0.5*( POP_STRUC.POPULATION(N).springK + POP_STRUC.POPULATION(N-1).springK );
        dx_1 = POP_STRUC.POPULATION(N+1).caCoordVector - POP_STRUC.POPULATION(N).caCoordVector;
        dx_2 = POP_STRUC.POPULATION(N).caCoordVector - POP_STRUC.POPULATION(N-1).caCoordVector;
        F_ela = tuo*( k1*norm(dx_1)-k2*norm(dx_2) );
        %F_ela = F_ela';
    end
    %F_pro
    %F_ela
    F_neb = (F_pro + F_ela) .*ORG_STRUC.whetherConstraint;
    
    %-----------------------------------------------%
    %------------ CI-image & DI-image  --------------%
    
    if (POP_STRUC.step >= ORG_STRUC.startCIDIStep) && (abs(ORG_STRUC.optMethodCIDI) > 0) && N>1 && N<numImages
        %
        % Automatic + Single CI-Image
        %
        if     isempty(ORG_STRUC.pickupImages) && ( ORG_STRUC.optMethodCIDI== 1 )
            if isCIDIImage(N)==max(isCIDIImage) && isCIDIImage(N)>0  % Only ONE image is selected!!
                disp(['==== Applying an automatic CI-NEB process at Image ' num2str(N) ' ! ...']);
                ORG_STRUC.whichCI=N;
                F_pro = ( F - ( 2*tuo*(tuo'*F) ) ) ;  %----- the Climbing-Image force
                F_ela(:) = 0;
            end
        end
        % Manual + Single CI-Image
        %ORG_STRUC.pickupImages
        if    ~isempty(ORG_STRUC.pickupImages) && ( ORG_STRUC.optMethodCIDI== 1 )
            if isCIDIImage(ORG_STRUC.pickupImages(1))<=0
                disp(' ');
                disp('ERROR: The set Image is not a Transition state, please check your input file! VCNEB Code quite with error ...');
                disp(' ');
                exit;
            elseif ORG_STRUC.pickupImages(1)==N
                disp(['==== Applying an manual CI-NEB process at Image ' num2str(N) ' ! ...']);
                ORG_STRUC.whichCI=N;
                F_pro = ( F - ( 2*tuo*(tuo'*F) ) ) ;  %----- the Climbing-Image force
                F_ela(:) = 0;
            end
        end
        %
        % Automatic + Single DI-Image
        %
        if     isempty(ORG_STRUC.pickupImages) && ( ORG_STRUC.optMethodCIDI==-1 )
            if isCIDIImage(N)==min(isCIDIImage) && isCIDIImage(N)<0
                disp(['==== Applying an automatic DI-NEB process at Image ' num2str(N) ' ! ...']);
                ORG_STRUC.whichDI=N;
                F_pro = F ;                           %----- the Downing-Image force
                F_ela(:) = 0;
            end
        end
        % Manual + Single DI-Image
        if    ~isempty(ORG_STRUC.pickupImages) && ( ORG_STRUC.optMethodCIDI==-1 )
            if isCIDIImage(ORG_STRUC.pickupImages(1))>=0
                disp(' ');
                disp('ERROR: The set Image is not a Local minimal, please check your input file! VCNEB Code quite with error ...');
                disp(' ');
                exit;
            elseif  (ORG_STRUC.pickupImages(1)==N)
                disp(['==== Applying an automatic DI-NEB process at Image ' num2str(N) ' ! ...']);
                ORG_STRUC.whichDI=N;
                F_pro = F ;                           %----- the Downing-Image force
                F_ela(:) = 0;
            end
        end
        
        
        %
        % Manual + Multi-CI/DI-Image
        %
        if  ~isempty(ORG_STRUC.pickupImages) && ( abs(ORG_STRUC.optMethodCIDI)>1 )
            
            disp('---Running multiCIDI-NEB process!--');
            for iC = 1:length(ORG_STRUC.pickupImages)
                isCIDIImage(ORG_STRUC.pickupImages(iC))
                if isCIDIImage(ORG_STRUC.pickupImages(iC))==0
                    disp(' ');
                    disp('ERROR: The set Image is not a Transition state or Local minimal, please check your input file! VCNEB Code quite with error ...');
                    disp(' ');
                    exit;
                elseif  ORG_STRUC.pickupImages(iC)==N
                    if         isCIDIImage(N)>0
                        disp(['==== Applying a CI of multiCIDI-NEB process at Image ' num2str(N) ' ! ...']);
                        if isempty(find(ORG_STRUC.whichCI==N))
                            ORG_STRUC.whichCI(end+1)=N;
                        end
                        
                        F_pro = ( F - ( 2*tuo*(tuo'*F) ) ) ;  %----- the Climbing-Image force
                        F_ela(:) = 0;
                    elseif     isCIDIImage(N)<0
                        disp(['==== Applying a DI of multiCIDI-NEB process at Image ' num2str(N) ' ! ...']);
                        if isempty(find(ORG_STRUC.whichDI==N))
                            ORG_STRUC.whichDI(end+1)=N;
                        end
                        F_pro = F ;                           %----- the Downing-Image force
                        F_ela(:) = 0;
                    end
                    
                end
            end
            
        end
    end
    %-----------------------------------------------%
    
    F_pro = F_pro .*ORG_STRUC.whetherConstraint;
    F_ela = F_ela .*ORG_STRUC.whetherConstraint;
    F_neb = F_pro + F_ela;
    
    POP_STRUC.POPULATION(N).F_neb = F_neb;
    POP_STRUC.POPULATION(N).F_pro = F_pro;
    POP_STRUC.POPULATION(N).F_ela = F_ela;
    
    %-----------------------------------------------------------------------------%
end



