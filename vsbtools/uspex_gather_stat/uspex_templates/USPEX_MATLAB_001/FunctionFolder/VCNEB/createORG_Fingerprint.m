function createORG_Fingerprint(inputFile)

global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fingerprints:
%[nothing, RmaxFing] = unix (['./getStuff ' inputFile ' RmaxFing 1']);
RmaxFing =  python_uspex(getPy, ['-f ' inputFile ' -b RmaxFing -c 1']);
if ~isempty(RmaxFing)
    ORG_STRUC.RmaxFing = str2num(RmaxFing);
end
if ORG_STRUC.RmaxFing > 0
    ORG_STRUC.doFing = 1;
else
    ORG_STRUC.doFing = 0;
end

%[nothing, deltaFing] = unix (['./getStuff ' inputFile ' deltaFing 1']);
deltaFing = python_uspex(getPy, ['-f ' inputFile ' -b deltaFing -c 1']);
if ~isempty(deltaFing)
    ORG_STRUC.deltaFing = str2num(deltaFing);
end

%[nothing, sigmaFing] = unix (['./getStuff ' inputFile ' sigmaFing 1']);
sigmaFing = python_uspex(getPy, ['-f ' inputFile ' -b sigmaFing -c 1']);
if ~isempty(sigmaFing)
    ORG_STRUC.sigmaFing = str2num(sigmaFing);
end

%[nothing, toleranceFing] = unix (['./getStuff ' inputFile ' toleranceFing 1']);
toleranceFing = python_uspex(getPy, ['-f ' inputFile ' -b toleranceFing -c 1']);
if ~isempty(toleranceFing)
    ORG_STRUC.toleranceFing = str2num(toleranceFing);
elseif ORG_STRUC.molecule
    ORG_STRUC.toleranceFing = 0.05;
end

% weight needed for normalisation of the cosine distance between fingerprints
if 1
    numIons = ORG_STRUC.numIons;
    L = length(numIons);
    S = 0;
    ORG_STRUC.weight = zeros(L*L,1);
    for i = 1:L
        for j = 1:L
            ind = (i-1)*L+j;
            ORG_STRUC.weight(ind) = (numIons(i)*numIons(j));
            S = S + (numIons(i)*numIons(j));
        end
    end
    
    ORG_STRUC.weight = ORG_STRUC.weight/S;
end
