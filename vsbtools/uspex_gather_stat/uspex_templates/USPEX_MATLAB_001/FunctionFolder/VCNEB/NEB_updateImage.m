function NEB_updateImage()


format long

global ORG_STRUC
global POP_STRUC

%--------     System Parameters      ------%
numImages = ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);
dimension  = ORG_STRUC.dimension;
numDimension = 3*(sumIons+dimension);

PressTarget = ORG_STRUC.ExternalPressure;
%-------------------------------------------------------------------------%

% Press and Volume
for i=1:numImages
    POP_STRUC.POPULATION(i).Press    = trace( POP_STRUC.POPULATION(i).cellStressMatirx/ORG_STRUC.KBar2GPa )/3;
    POP_STRUC.POPULATION(i).cellStressMatirx/ORG_STRUC.KBar2GPa;
    POP_STRUC.POPULATION(i).VOLUME   = abs( det( POP_STRUC.POPULATION(i).LATTICE ) );
    POP_STRUC.POPULATION(i).Energy   = POP_STRUC.POPULATION(i).Enthalpy - PressTarget*POP_STRUC.POPULATION(i).VOLUME/160.217733;
    %  eV = eV + GPa*A^3/160.217733
end

for i=1:numImages
    L = POP_STRUC.POPULATION(i).LATTICE;
    %  Force on the Atoms
    POP_STRUC.POPULATION(i).aCaForceMatrix = POP_STRUC.POPULATION(i).atomForcesMatrix ;                % eV/A
    %    disp('aCaForceMatrix : ')
    %    POP_STRUC.POPULATION(i).aCaForceMatrix
    POP_STRUC.POPULATION(i).aFrForceMatrix = ( (L'*L)*L^(-1)*(POP_STRUC.POPULATION(i).aCaForceMatrix') )';       % eV
    %  Stress on the Cell
    %  GPa*A^3   -->  eV
    % if 1==vcnebopt.PressType
    PressTensor=eye(3)*PressTarget;
    %  else
    %    PressTensor=PressTarget;
    %  end
    POP_STRUC.POPULATION(i).cFrStressMatrix =  POP_STRUC.POPULATION(i).VOLUME*( POP_STRUC.POPULATION(i).cellStressMatirx /ORG_STRUC.KBar2GPa - PressTensor )/(160.217733);
    tmp =( ( L^(-1) )'*POP_STRUC.POPULATION(i).cFrStressMatrix ) ;
    POP_STRUC.POPULATION(i).cCaStressMatrix = tmp - POP_STRUC.POPULATION(i).COORDINATES'*POP_STRUC.POPULATION(i).aCaForceMatrix;
    
    %------------ Setting X_image & u_image --------------%
    %---[ Cell lattice : 9  + Atom Postion : 3*sumIons  ]
    %---[  1  2  3  ]
    %---|  4  5  6  ]
    %---[  7  8  9  ]
    % For the Cartesical Coordination
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-2)-2) = L(1,1);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-2)-1) = L(1,2);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-2)-0) = L(1,3);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-1)-2) = L(2,1);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-1)-1) = L(2,2);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-1)-0) = L(2,3);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-0)-2) = L(3,1);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-0)-1) = L(3,2);
    POP_STRUC.POPULATION(i).caCoordVector(3*(dimension-0)-0) = L(3,3);
    % For the Fractional Coordination
    U = POP_STRUC.POPULATION(i).caCoordVector(1:3*dimension);
    POP_STRUC.POPULATION(i).frCoordVector(1:3) = U(1:3)/norm( U(1:3));
    POP_STRUC.POPULATION(i).frCoordVector(4:6) = U(4:6)/norm( U(4:6));
    POP_STRUC.POPULATION(i).frCoordVector(7:9) = U(7:9)/norm( U(7:9));
    
    for j = 1+dimension:sumIons+dimension    % For the Cartesical Coordination
        sx = POP_STRUC.POPULATION(i).CARTECOORDS(j-dimension,1);   %  bohr --> A
        sy = POP_STRUC.POPULATION(i).CARTECOORDS(j-dimension,2);   %  bohr --> A
        sz = POP_STRUC.POPULATION(i).CARTECOORDS(j-dimension,3);   %  bohr --> A
        POP_STRUC.POPULATION(i).caCoordVector( (j-1)*3+1:j*3) = [sx, sy, sz];
    end
    for j = 1+dimension:sumIons+dimension    % For the Fractional Coordination
        sx = POP_STRUC.POPULATION(i).COORDINATES(j-dimension,1);   %  bohr --> A
        sy = POP_STRUC.POPULATION(i).COORDINATES(j-dimension,2);   %  bohr --> A
        sz = POP_STRUC.POPULATION(i).COORDINATES(j-dimension,3);   %  bohr --> A
        POP_STRUC.POPULATION(i).frCoordVector( (j-1)*3+1:j*3) = [sx, sy, sz];
    end
end

%--------------------------------------------------------------------
%Force = [Cell Force*9 + Natom Force * 3sumIons ]
%---[  1  2  3  ]
%---|  4  5  6  ]
%---[  7  8  9  ]
for i = 1:numImages
    for j = 1+dimension:sumIons+dimension
        POP_STRUC.POPULATION(i).totFrForceVector( (3*(j-1)+1):(3*j)) = POP_STRUC.POPULATION(i).aFrForceMatrix(j-dimension,:);
        POP_STRUC.POPULATION(i).totCaForceVector( (3*(j-1)+1):(3*j)) = POP_STRUC.POPULATION(i).aCaForceMatrix(j-dimension,:);
    end
    j = dimension;
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-2)-2) = POP_STRUC.POPULATION(i).cFrStressMatrix(1,1);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-2)-1) = POP_STRUC.POPULATION(i).cFrStressMatrix(1,2);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-2)-0) = POP_STRUC.POPULATION(i).cFrStressMatrix(1,3);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-1)-2) = POP_STRUC.POPULATION(i).cFrStressMatrix(2,1);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-1)-1) = POP_STRUC.POPULATION(i).cFrStressMatrix(2,2);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-1)-0) = POP_STRUC.POPULATION(i).cFrStressMatrix(2,3);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-0)-2) = POP_STRUC.POPULATION(i).cFrStressMatrix(3,1);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-0)-1) = POP_STRUC.POPULATION(i).cFrStressMatrix(3,2);
    POP_STRUC.POPULATION(i).totFrForceVector( 3*(j-0)-0) = POP_STRUC.POPULATION(i).cFrStressMatrix(3,3);
    
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-2)-2) = POP_STRUC.POPULATION(i).cCaStressMatrix(1,1);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-2)-1) = POP_STRUC.POPULATION(i).cCaStressMatrix(1,2);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-2)-0) = POP_STRUC.POPULATION(i).cCaStressMatrix(1,3);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-1)-2) = POP_STRUC.POPULATION(i).cCaStressMatrix(2,1);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-1)-1) = POP_STRUC.POPULATION(i).cCaStressMatrix(2,2);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-1)-0) = POP_STRUC.POPULATION(i).cCaStressMatrix(2,3);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-0)-2) = POP_STRUC.POPULATION(i).cCaStressMatrix(3,1);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-0)-1) = POP_STRUC.POPULATION(i).cCaStressMatrix(3,2);
    POP_STRUC.POPULATION(i).totCaForceVector( 3*(j-0)-0) = POP_STRUC.POPULATION(i).cCaStressMatrix(3,3);
    %       disp('totCaForceVector : ')
    %       POP_STRUC.POPULATION(i).totCaForceVector
end
%------------------------------

%-------  Display the Image INfo   --------------%
%for i = 1:numImages

%end
