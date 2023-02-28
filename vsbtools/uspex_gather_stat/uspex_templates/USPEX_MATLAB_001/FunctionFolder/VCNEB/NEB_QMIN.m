function  NEB_QMIN()

global ORG_STRUC
global POP_STRUC

%--------     System Parameters      ------%
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);
dimension  = ORG_STRUC.dimension;
numDimension=3*(sumIons+dimension);

dt=ORG_STRUC.dt;

%-- Strucutre Optimization Algorithm  ---
%    1: Quick-MIN / QM
%    2: FIRE                 (default)
%    3: BFGS                      (not yet)
%    4: Fast Method for Cell Opt  (not yet)

if ORG_STRUC.optimizerType >2
    disp(['ERROR: Optimize algorithm other than Quick-MIN and FIRE is not implemented yet!']);
    exit;
end



for N = 1:numImages
    alpha = 1;
    move_rangeMax = 0;
    g_pro = POP_STRUC.POPULATION(N).F_pro;
    g_ela = POP_STRUC.POPULATION(N).F_ela;
    
    %-----  Steepest-Descent Algorithm
    if 1==ORG_STRUC.optimizerType   % more need to accelary !!!!
        if N==1 || N==numImages
            dx_QMIN = NEB_algoFire(N);
        else
            NEB_algoFire(N);
            dx_pro = 1/2*g_pro*dt*dt;
            dx_ela = 1/2*g_ela*dt*dt;
            dx_QMIN=dx_pro+dx_ela;
        end
        %dx_QMIN = NEB_trustRadius(dx_pro) + dx_ela;
    end
    %----- Fire Algorithm (PRL)
    if 2==ORG_STRUC.optimizerType
        dx_QMIN = NEB_algoFire(N);
    end
    
    %-----    Conjugate Gradients Alogrithm
    if 3==ORG_STRUC.optimizerType
        % Not implemented Well!!~~
        if 1==POP_STRUC.step
            [H_phon, alpha] = NEB_initHessian( H*g, 1, move_rangeMax, N);
            dx_QMIN = -1*H_phon*(-g);
            Hessian(:,:,N) = H_phon;
        end
        if POP_STRUC.step>1
            [d, H_phon] = NEB_cgvar(N);
            dx_QMIN = - 1*H_phon*(-d) ;
        end
    end
    %-----     BFGS Algorithm
    if 4==ORG_STRUC.optimizerType
        %-------------------------------------------------%
        if 1==POP_STRUC.step
            [H_phon, alpha] = NEB_initHessian( H*g, 1, move_rangeMax, N);
            dx_QMIN = -1*H_phon*(-g);
            Hessian(:,:,N) = H_phon;
        end
        if POP_STRUC.step>1
            H_phon = NEB_bfgsvar(N);
            dx_QMIN = -1*H_phon*(-g);
            Hessian(:,:,N) = H_phon;
            display(dx_QMIN);
        end
    end
    dX = NEB_trustRadius(dx_QMIN, N);
    POP_STRUC.POPULATION(N).X_move = dX;
end


