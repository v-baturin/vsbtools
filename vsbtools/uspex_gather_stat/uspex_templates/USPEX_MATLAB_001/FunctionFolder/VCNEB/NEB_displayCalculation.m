function NEB_displayCalculation()



%-----------------------------------------------%
global ORG_STRUC
global POP_STRUC


numImages = ORG_STRUC.numImages;
numIons = ORG_STRUC.numIons;
sumIons = sum(ORG_STRUC.numIons);

NEB_findSymmetry();

cd(ORG_STRUC.homePath);
fpath = [ORG_STRUC.resFolder '/' ORG_STRUC.log_file];
fp = fopen(fpath, 'a');
%-----------------------------------------------%
if 1==ORG_STRUC.CalcType
    fprintf(fp, '=     VCNEB Calculation Detail Information :     STEP %6d        \n=\n=', POP_STRUC.generation);
else
    fprintf(fp, '=   Relaxation Calculation Detail Information :     Step: %6d \n=\n=', POP_STRUC.generation);
end

%-- Calculation Details
fprintf(fp, '         Number of images  =    %2d \n', numImages);

if   0==mod(POP_STRUC.step, ORG_STRUC.PrintStep) || (POP_STRUC.step==ORG_STRUC.numSteps)
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
        str = '.GULP.';
    end
    fprintf(fp, '   DFT Software Interface  =  %s \n',  str);
    %----------------------------------------------------------------------
    
    fprintf(fp, '        Extenral Pressure  =   %5.2f GPa \n',ORG_STRUC.ExternalPressure  );
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
        fprintf(fp, '     Variable Path Length  = %6.3f A\n',  ORG_STRUC.VarPathLength);
    else
        fprintf(fp, '           Variable Image  =  .False. \n');
    end
    %----------------------------------------------------------------------
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
    if  1==ORG_STRUC.optFreezing
        str = '.True.';
    else
        str = '.False.';
    end
    fprintf(fp, '                 Freezing  =  %s \n',  str);
    %----------------------------------------------------------------------
    fprintf(fp, '  Force RMS Halt Criteria  =   %5.3f eV/A\n',  ORG_STRUC.ConvThreshold);
    
    
    if abs(ORG_STRUC.optMethodCIDI) > 0 && 1==ORG_STRUC.CalcType
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
    %%%%%%%%---------------------------------------------------------------------------------------
    for i = 1:numImages
        % Lattice Parameter
        fprintf(fp,'\n                           ==================================================\n');
        if POP_STRUC.POPULATION(i).freezing==1
            fprintf(fp,'                           =    Image Number: %3d,    ( Freezing=.TRUE. )   =', i);
        else
            fprintf(fp,'                           =    Image Number: %3d,    ( Freezing=.FALSE. )  =', i);
        end
        fprintf(fp,'\n                           ==================================================\n');
        fprintf(fp, 'Space group: %s \n', POP_STRUC.POPULATION(i).spacegroupSymbol);
        
        %--------
        fprintf(fp, 'Lattice Parameters:\n');
        Lattice = latConverter(POP_STRUC.POPULATION(i).LATTICE);
        if size(Lattice,1) == 1
            Lattice = Lattice';
        end
        Lattice(4:6) = Lattice(4:6)*180/pi;
        fprintf(fp, '        a = %7.4f,    b = %7.4f,     c = %7.4f    (A)\n', Lattice(1:3));
        fprintf(fp, '    alpha = %7.3f, beta = %7.3f, gamma = %7.3f    (degree)\n', Lattice(4:6));
        fprintf(fp, '\n\n');
        
        
        fprintf(fp, 'External Pressure (kBar) :  \n');
        fprintf(fp, '   %12.6f\n\n', trace(POP_STRUC.POPULATION(i).cellStressMatirx)/3);
        fprintf(fp, '          Cell Vector Matrix (A),               Force on cell     (kBar) ,             Lattice Strain  (A):\n');
        for iC = 1:3
            fprintf(fp, '       %10.6f %10.6f %10.6f,  ', POP_STRUC.POPULATION(i).LATTICE(iC,:));
            fprintf(fp, '  %10.4f %10.4f %10.4f,  ', POP_STRUC.POPULATION(i).cellStressMatirx(iC,:));
            fprintf(fp, '  %10.6f %10.6f %10.6f\n', POP_STRUC.POPULATION(i).LATTICE_move(iC,:));
        end
        
        
        
        % Atomic Information
        fprintf(fp,'\n\n       Atomic Fractional Coordinates,      Forces acting on Atoms (eV/A*Vol^1/3)      Displacement of Atoms (A)\n');
        coordLoop = 1;
        for iC = 1 : length(numIons)
            for j = 1 : numIons(iC)
                fprintf(fp, '%4s   %10.6f %10.6f %10.6f,  ', megaDoof(ORG_STRUC.atomType(iC)), POP_STRUC.POPULATION(i).COORDINATES(coordLoop,:));
                fprintf(fp, '  %10.6f %10.6f %10.6f,  ', POP_STRUC.POPULATION(i).atomForcesMatrix(coordLoop,:));
                fprintf(fp, '  %10.6f %10.6f %10.6f\n', POP_STRUC.POPULATION(i).COORDINATES_move(coordLoop,:));
                
                coordLoop = coordLoop + 1;
            end
        end
        
        
        if 1==ORG_STRUC.CalcType
            
            cellFneb= reshape( POP_STRUC.POPULATION(i).F_neb(1:9), 3,3 );
            cellFpro= reshape( POP_STRUC.POPULATION(i).F_pro(1:9), 3,3 );
            cellFela= reshape( POP_STRUC.POPULATION(i).F_ela(1:9), 3,3 );
            
            %fprintf(fp,'\n\nDetial of VCNEB forces ( eV/A) :\n');
            fprintf(fp, ['\n\n' alignLine('VCNEB forces on Lattice: ( eV/A*Vol^1/3 )') '\n']);
            
            fprintf(fp,'              VCENB Total Forces,              VCNEB Project Forces,             VCNEB Elastic Forces\n');
            for iC = 1:3
                fprintf(fp, '     [%10.6f %10.6f %10.6f], [%10.6f %10.6f %10.6f], [%10.6f %10.6f %10.6f]\n',...
                    cellFneb(iC,:),cellFpro(iC,:),cellFela(iC,:));
            end
            fprintf(fp,'\nLattice VCNEB forces Maximum & RootMeanSquare (eV/A)\n' );
            fprintf(fp,'    VCNEB Forces Max:       %10.6f               %10.6f               %10.6f\n',...
                POP_STRUC.POPULATION(i).errcFnebMax,POP_STRUC.POPULATION(i).errcFproMax,POP_STRUC.POPULATION(i).errcFelaMax);
            fprintf(fp,'    VCNEB Forces RMS:       %10.6f               %10.6f               %10.6f\n',...
                POP_STRUC.POPULATION(i).errcFnebRms,POP_STRUC.POPULATION(i).errcFproRms,POP_STRUC.POPULATION(i).errcFelaRms);
            
            
            atomFneb= reshape( POP_STRUC.POPULATION(i).F_neb(10:end), sumIons,3 );
            atomFpro= reshape( POP_STRUC.POPULATION(i).F_pro(10:end), sumIons,3 );
            atomFela= reshape( POP_STRUC.POPULATION(i).F_ela(10:end), sumIons,3 );
            %VCNEB forces on atoms
            fprintf(fp, ['\n' alignLine('VCNEB forces on Atoms: ( eV/A*Vol^1/3)') '\n']);
            fprintf(fp,'              VCENB Total Forces,              VCNEB Project Forces,             VCNEB Elastic Forces\n');
            
            coordLoop = 1;
            for iC = 1 : length(numIons)
                for j = 1 : numIons(iC)
                    fprintf(fp, '%4s [%10.6f %10.6f %10.6f], [%10.6f %10.6f %10.6f], [%10.6f %10.6f %10.6f]\n',...
                        megaDoof(ORG_STRUC.atomType(iC)),atomFneb(coordLoop,:),atomFpro(coordLoop,:),atomFela(coordLoop,:));
                    coordLoop = coordLoop + 1;
                end
            end
            fprintf(fp,'\nAtomic VCNEB Forces Maximum & RootMeanSquare (eV/A)\n' );
            fprintf(fp,'    VCNEB Forces Max:       %10.6f               %10.6f               %10.6f\n',...
                POP_STRUC.POPULATION(i).erraFnebMax,POP_STRUC.POPULATION(i).erraFproMax,POP_STRUC.POPULATION(i).erraFelaMax);
            fprintf(fp,'    VCNEB Forces RMS:       %10.6f               %10.6f               %10.6f\n',...
                POP_STRUC.POPULATION(i).erraFnebRms,POP_STRUC.POPULATION(i).erraFproRms,POP_STRUC.POPULATION(i).erraFelaRms);
            
        end
        % Fire Algorithm Info
        fprintf(fp, 'Optimize Algorithm Parameters (FIRE) : \n');
        fprintf(fp, '      Ncount = %3d,  dt = %8.6f,  a= %8.4f,  P= %8.6f \n',...
            POP_STRUC.POPULATION(i).Fire.Ncount, POP_STRUC.POPULATION(i).Fire.dt, POP_STRUC.POPULATION(i).Fire.a, POP_STRUC.POPULATION(i).Fire.P);
        fprintf(fp, 'FIRE velocity:\n');
        fprintf(fp, '%11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f \n', POP_STRUC.POPULATION(i).Fire.v);
        fprintf(fp, '\n');
        % JOB infomation
        
    end
    
end



fclose(fp);
