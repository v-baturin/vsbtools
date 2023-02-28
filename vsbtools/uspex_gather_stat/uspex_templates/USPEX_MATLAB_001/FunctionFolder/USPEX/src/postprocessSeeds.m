function goodPOSCARs = postprocessSeeds(POSCARs, sysAtomType, sysNumIons, STDMOL)


goodPOSCARs = '';

numOfGood = 0;
for ID = 1:length(POSCARs)
    scale_Seeds    = POSCARs(ID).scale;
    lattice_Seeds  = POSCARs(ID).lattice;
    atomType_Seeds = POSCARs(ID).atomType;
    numMols_Seeds  = POSCARs(ID).numMols;
    numIons_Seeds  = POSCARs(ID).numIons;
    coord_Seeds    = POSCARs(ID).coord;
    atomType       = GetElement(length(numIons_Seeds), atomType_Seeds);   

    %-- check the consistent of AtomType
    isAtomTypeOK=1;
     
    if length(numIons_Seeds) ~= length(atomType) %more atomTypes
        USPEXmessage(552, '', 0);
        isAtomTypeOK = 0;
    elseif sum(ismember(atomType, sysAtomType))<length(numIons_Seeds)
        USPEXmessage(553, '', 0);
        isAtomTypeOK = 0;
    end

    if  isAtomTypeOK == 0
        warningStr = ['Seeds : Element types are inconsistent in Seeds-' num2str(ID) ];
        USPEXmessage(0, warningStr, 0);

        POSCARs(ID)

        continue;
    end

    %--- fill all the missing atomType
    numIons = zeros(1, length(sysAtomType));
    for i=1:length(sysAtomType)
        for j=1:length(numIons_Seeds)
            if atomType(j)==sysAtomType(i)
                numIons(i) = numIons_Seeds(j);
            end
        end
    end

    %--- maping the atomic positions
    
    coord = zeros(size(coord_Seeds));
    startIons=0;
    for i = 1:length(sysAtomType)
        for j = 1:length(atomType)
            if sysAtomType(i)==atomType(j)
                if j==1
                    coord(startIons+1:startIons+numIons_Seeds(j),:)= ...
                                       coord_Seeds(1:numIons_Seeds(1),:);
                    startIons = startIons+numIons_Seeds(j);
                else
                    coord(startIons+1:startIons+numIons_Seeds(j),:)= ...
                    coord_Seeds( sum(numIons_Seeds(1:j-1))+1:sum(numIons_Seeds(1:j)),: );
                    startIons = startIons+numIons_Seeds(j);
                end
            end
        end
    end
    coord = coord*scale_Seeds;

    % Mole
    numMols    = 1;   % fix the bug,  mean nothing here.
    numBlocks = numIons/sysNumIons;
    if ~isempty(STDMOL)
        comp=zeros(size(sysNumIons,2), length(sysAtomType));
        for i=1:size(sysNumIons,2)
            for j=1:length(STDMOL(i).types)
                for k=1:length(sysAtomType)
                    if STDMOL(i).types(j)==k
                        comp(i,k)=comp(i,k)+1;
                    end
                end
            end
        end
        if ~isempty(numMols_Seeds) ...
           && length(numMols_Seeds) ==size( comp, 1 ) ...
           && norm(numIons - numMols_Seeds*comp ) < 1e-3
            numMols = numMols_Seeds;
        else
            numMols = round(numIons/comp);
        end
        numBlocks = numMols/sysNumIons;
        
    end

    numOfGood = numOfGood + 1;

    goodPOSCARs(numOfGood).scale    = 1.0;
    goodPOSCARs(numOfGood).atomType = atomType;
    goodPOSCARs(numOfGood).numIons  = numIons;
    goodPOSCARs(numOfGood).lattice  = lattice_Seeds;
    goodPOSCARs(numOfGood).coord    = coord;
    goodPOSCARs(numOfGood).numMols  = numMols;
    goodPOSCARs(numOfGood).numBlocks= numBlocks;

end


