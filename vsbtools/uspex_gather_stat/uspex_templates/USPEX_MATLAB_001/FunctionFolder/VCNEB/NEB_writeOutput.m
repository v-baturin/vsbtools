function NEB_writeOutput(step)

global ORG_STRUC
global POP_STRUC


%--------     System Parameters      ------%
numImages = ORG_STRUC.numImages;
sumIons = sum(ORG_STRUC.numIons);

%-----------------------------------------------%
cd(ORG_STRUC.homePath);
fpath = [ORG_STRUC.resFolder '/VCNEBReports'];

if -1==step
    fp = fopen(fpath, 'a');
    fprintf(fp, [alignLine('VCNEB Method Calculation Brief Report') '\n']);
    fprintf(fp, '\n');
    %fprintf(fp, '         Number of images  =    %2d \n', numImages);
    %-----------------------------------------------------------------------
    if     1==ORG_STRUC.CalcType
        fprintf(fp, '         Calculation Type  =  .VCNEB.  \n');
    elseif 2==ORG_STRUC.CalcType
        fprintf(fp, '         Calculation Type  =  .Relaxation. \n');
    end
    %----------------------------------------------------------------------
    if       1==ORG_STRUC.abinitioCode
        str = '.VASP.';
    elseif   8==ORG_STRUC.abinitioCode
        str = '.Quantum Espresso.';
    elseif   3==ORG_STRUC.abinitioCode
        str = '.ABINIT.';
    end
    fprintf(fp, '   DFT Software Interface  =  %s \n',  str);
    %----------------------------------------------------------------------
    
    fprintf(fp, '         Extenral Pressure =   %5.2f GPa \n',ORG_STRUC.ExternalPressure  );
    %----------------------------------------------------------------------
    if       1==ORG_STRUC.optimizerType
        str = '.Steep Descent.';
    elseif   2==ORG_STRUC.optimizerType
        str = '.Fire.';
    elseif   3==ORG_STRUC.optimizerType
        str = '.Verlet.';
    elseif   2==ORG_STRUC.optimizerType
        str = '.BFGS.';
    end
    fprintf(fp, '   Optimisation Algorithm  =  %s \n',  str);
    %----------------------------------------------------------------------
    %   if   isempty(ORG_STRUC.pickupImages)
    %       str = '.True.';
    %       fprintf(fp, '  No Peak Images specified in this calculation \n');
    %   else
    %       fprintf(fp, ' Peak Image = %s \n', num2str(ORG_STRUC.pickupImages));
    %       fprintf(fp, ' Method (1: Climbing; 2: Downing) \n');
    %       fprintf(fp, ' Method     = %s \n', num2str(ORG_STRUC.Method));
    %   end
    %----------------------------------------------------------------------
    if  2==ORG_STRUC.optRelaxType
        str = '.True.';
    else
        str = '.False.';
    end
    fprintf(fp, '  Lattice Cell Relaxation  =  %s \n',  str);
    
    %----------------------------------------------------------------------
    
    if 1==ORG_STRUC.CalcType
        if 1==ORG_STRUC.optVarK
            fprintf(fp, ' Variable Spring Constant  =  .True. \n');
            fprintf(fp, '                    --  k-min  =  %3.1f eV/A\n',  ORG_STRUC.K_min);
            fprintf(fp, '                    --  k-max  =  %3.1f eV/A\n',  ORG_STRUC.K_max);
        else
            fprintf(fp, ' Variable Spring Constant  =  .False. \n');
            fprintf(fp, '               --  k-constant  =   %5.1f eV/A\n',  ORG_STRUC.Kconstant);
        end
    end
    
    %----------------------------------------------------------------------
    if 1==ORG_STRUC.optVarImage
        fprintf(fp, '    Variable Image Number  =  .True. \n');
        fprintf(fp, '    Variable Path Length   = %6.3f A\n',  ORG_STRUC.VarPathLength);
    else
        fprintf(fp, '           Variable Image  =  .False. \n');
    end
    %----------------------------------------------------------------------
    
    if ORG_STRUC.optMethodCIDI~=0
        if ORG_STRUC.optMethodCIDI==1
            fprintf(fp, '  Climbing/Downing Images  =  .CI.\n');
        elseif ORG_STRUC.optMethodCIDI==-1
            fprintf(fp, '  Climbing/Downing Images  =  .DI.\n');
        elseif ORG_STRUC.optMethodCIDI>1
            fprintf(fp, '  Climbing/Downing Images  =  .multi-CI/DI.\n');
        end
        if  ~isempty(ORG_STRUC.pickupImages)
            fprintf(fp, '  picked CI/DI Images No.  =   ');
            fprintf(fp, '%2d, ', ORG_STRUC.pickupImages);
            fprintf(fp, '\n');
        end
    end
    %----
    if  1==ORG_STRUC.optFreezing
        str = '.True.';
    else
        str = '.False.';
    end
    fprintf(fp, '                 Freezing  =  %s \n',  str);
    
    fprintf(fp, '  Force RMS Halt Criteria  =   %5.3f eV/A\n',  ORG_STRUC.ConvThreshold);
    
    fprintf(fp, '\n');
    fclose(fp);
    
    if ORG_STRUC.numSteps~=-1
        return;
    end
end

isTS=NEB_isTransitionState();
if (POP_STRUC.step>=0)
    
    fp = fopen(fpath, 'a');
    if POP_STRUC.step==0
        if 2==ORG_STRUC.CalcType
            fprintf(fp, [alignLine('Relaxation step : initial') '\n']);
        else
            fprintf(fp, [alignLine('VCNEB step : initial') '\n']);
        end
    else
        if 2==ORG_STRUC.CalcType
            fprintf(fp, [alignLine( sprintf('Relaxation step : %5d', POP_STRUC.step) ) '\n'] );
        else
            fprintf(fp, [alignLine( sprintf('VCNEB step :  %5d', POP_STRUC.step) ) '\n'] );
        end
    end
    fprintf(fp, '                      Cell Forces RMS(eV/A)      Atom Forces RMS(eV/A)              \n');
    fprintf(fp, 'Image Energy-(eV)    Total, Project, Elastic    Total, Project, Elastic  F   Dist   SpaceG\n');
    for i = 1:numImages
        H(i)=POP_STRUC.POPULATION(i).Enthalpy;
        erraFnebRms(i)=POP_STRUC.POPULATION(i).erraFnebRms;
        erraFproRms(i)=POP_STRUC.POPULATION(i).erraFproRms;
        erraFelaRms(i)=POP_STRUC.POPULATION(i).erraFelaRms;
        errcFnebRms(i)=POP_STRUC.POPULATION(i).errcFnebRms;
        errcFproRms(i)=POP_STRUC.POPULATION(i).errcFproRms;
        errcFelaRms(i)=POP_STRUC.POPULATION(i).errcFelaRms;
        freezing(i)=POP_STRUC.POPULATION(i).freezing;
        L_path(i)=POP_STRUC.POPULATION(i).pathLength;
        spacegroupNumber(i)=POP_STRUC.POPULATION(i).spacegroupNumber;
        fprintf(fp, '%3d    %10.4f  ', i,  H(i));
        fprintf(fp, '[%7.3f,%7.3f,%7.3f]  [%7.3f,%7.3f,%7.3f]  %1d  %6.4f  %4d', ...
            errcFnebRms(i),  errcFproRms(i),   errcFelaRms(i),  erraFnebRms(i),  erraFproRms(i),   erraFelaRms(i), freezing(i), L_path(i), spacegroupNumber(i) );
        if ORG_STRUC.CalcType==1
            if isTS(i)>0
                fprintf(fp,'    <-TS->\n');
            elseif isTS(i)<0
                fprintf(fp,'    ->LM<-\n');
            else
                fprintf(fp,'          \n');
            end
        else
            fprintf(fp,'          \n');
        end
    end
    
    if 1==ORG_STRUC.CalcType
        %-----------------------------------------%
        fprintf(fp, '\n   Activation energy (->) = %10.3f eV',   ...
            (max(H(2:numImages-1))-H(1)) );
        fprintf(fp, '\n   Activation energy (<-) = %10.3f eV\n', ...
            (max(H(2:numImages-1))-H(numImages)) );
        %--------------%
    end
    if abs(ORG_STRUC.optMethodCIDI) > 0 && (1==ORG_STRUC.CalcType)
        if       ORG_STRUC.optMethodCIDI == 1
            fprintf(fp, '\n   Climbing Image No. for Transition State :  %2d  \n', ORG_STRUC.whichCI);
        elseif   ORG_STRUC.optMethodCIDI ==-1
            fprintf(fp, '\n   Downing  Image No. for Metastable State :  %2d  \n', ORG_STRUC.whichDI);
        elseif   ORG_STRUC.optMethodCIDI > 1
            fprintf(fp, '\n   Multi-Climbing Image No. for Transition States :  ');
            fprintf(fp, '%2d, ', ORG_STRUC.whichCI);
            
            fprintf(fp, '\n   Multi-Downing  Image No. for Metastable States :  ');
            fprintf(fp, '%2d, ', ORG_STRUC.whichDI); fprintf(fp, '\n');
        end
        
    end
    
    fclose(fp);
    
    %--------------------------------------------------------------------
    if ORG_STRUC.numSteps~=-1
        fpath = [POP_STRUC.resFolder '/Energy'];
        fp = fopen(fpath, 'w');
        if 1==ORG_STRUC.CalcType
            fprintf(fp,'VCNEB Calcualtion Step %4d\n\n', POP_STRUC.step);
        else
            fprintf(fp,'Relaxation Calcualtion Step %4d\n\n', POP_STRUC.step);
        end
        fprintf(fp, 'Image  Energy-(eV)  Dist \n');
        for i = 1:numImages
            fprintf( fp, '%3d   %10.4f  %6.4f\n', i,  H(i), L_path(i) );
        end
        fclose(fp);
    end
end
%---------------------------------------------------------------------------------------------%
