function NEB_calcError()

%--------------------------------------%
global ORG_STRUC
global POP_STRUC

%--------------------------------------%
numImages =     ORG_STRUC.numImages;
sumIons  = sum(ORG_STRUC.numIons);

dimension  = ORG_STRUC.dimension;
numDimension = 3*(sumIons+dimension);


cellFnebMatrix=zeros(3,dimension);                          %vcneb forces on cell Lattice in Matrix form
cellFelaMatrix=zeros(3,dimension);
cellFproMatrix=zeros(3,dimension);

atomFnebMatrix=zeros(3,sumIons);                  %vcneb forces on atoms in Matrix form
atomFelaMatrix=zeros(3,sumIons);
atomFproMatrix=zeros(3,sumIons);
%--------------------------------------^
for N = 1:numImages
    F_neb=POP_STRUC.POPULATION(N).F_neb;
    F_pro=POP_STRUC.POPULATION(N).F_pro;
    F_ela=POP_STRUC.POPULATION(N).F_ela;
    %------------------------------------------
    cellFnebMatrix=reshape(F_neb(1:3*dimension),3,dimension)';  %  [ 1 2 3]
    cellFelaMatrix=reshape(F_ela(1:3*dimension),3,dimension)';  %  [ 4 5 6]
    cellFproMatrix=reshape(F_pro(1:3*dimension),3,dimension)';  %  [ 7 8 9]
    atomFnebMatrix=reshape(F_neb(3*dimension+1:end),3,sumIons); % [ 1 4 7 10 ...
    atomFelaMatrix=reshape(F_ela(3*dimension+1:end),3,sumIons); % [ 2 5 8 11 ...
    atomFproMatrix=reshape(F_pro(3*dimension+1:end),3,sumIons); % [ 3 6 9 12 ...
    
    
    L = POP_STRUC.POPULATION(N).LATTICE;
    factor=(L'*L)*L^(-1);
    
    errorCellNebMatrix=(factor^(-1)*(cellFnebMatrix));
    errorCellElaMatrix=(factor^(-1)*(cellFelaMatrix));
    errorCellProMatrix=(factor^(-1)*(cellFproMatrix));
    errorAtomNebMatrix=(factor^(-1)*(atomFnebMatrix));
    errorAtomElaMatrix=(factor^(-1)*(atomFelaMatrix));
    errorAtomProMatrix=(factor^(-1)*(atomFproMatrix));
    for j = 1:sumIons
        errorAtomNebRMS(j)=norm(errorAtomNebMatrix(:,j))/sqrt(3);
        errorAtomElaRMS(j)=norm(errorAtomElaMatrix(:,j))/sqrt(3);
        errorAtomProRMS(j)=norm(errorAtomProMatrix(:,j))/sqrt(3);
    end
    %%%  -- max Forces on cell Lattice and atoms
    POP_STRUC.POPULATION(N).errcFnebMax=max(max(errorCellNebMatrix ));
    POP_STRUC.POPULATION(N).errcFproMax=max(max(errorCellProMatrix ));
    POP_STRUC.POPULATION(N).errcFelaMax=max(max(errorCellElaMatrix ));
    POP_STRUC.POPULATION(N).erraFnebMax=max(errorAtomNebRMS);
    POP_STRUC.POPULATION(N).erraFproMax=max(errorAtomProRMS);
    POP_STRUC.POPULATION(N).erraFelaMax=max(errorAtomElaRMS);
    
    %%% RMS on cell Lattice and atoms
    POP_STRUC.POPULATION(N).errcFnebRms=norm(errorCellNebMatrix)/sqrt(3*dimension);
    POP_STRUC.POPULATION(N).errcFelaRms=norm(errorCellElaMatrix)/sqrt(3*dimension);
    POP_STRUC.POPULATION(N).errcFproRms=norm(errorCellProMatrix)/sqrt(3*dimension);
    POP_STRUC.POPULATION(N).erraFnebRms=norm(errorAtomNebMatrix)/sqrt(3*sumIons);
    POP_STRUC.POPULATION(N).erraFelaRms=norm(errorAtomElaMatrix)/sqrt(3*sumIons);
    POP_STRUC.POPULATION(N).erraFproRms=norm(errorAtomProMatrix)/sqrt(3*sumIons);
end

