function  POSCARs = splitSeedsContext(seedContent)


%The first line of seedContent is a comment

POSCARs='';

numOfPOSCAR = 0;

startline = 1;
while startline <= length( seedContent )-7
    startline = startline+1;
    
    numMols = [];
    scale   = [];
    atomType= '';
    numIons = [];
    lattice = [];
    coord = [];
    
    try
        numMols = numMolReadSeeds( char(seedContent(startline)) );
        % Read the MOL user setting.
        scale   = str2num( char(seedContent(startline + 1)) );
        if length(scale) ~= 1
            seedProcessError(numOfPOSCAR, startline + 1);
            continue;
        end
        
        lattice = str2num(char(seedContent(startline + 2:startline + 4)));
        if isempty(lattice)
            seedProcessError(numOfPOSCAR, startline + 2);
            continue;
        end
        
        atomType= ( char(seedContent(startline + 5)) );
        numIons = str2num( char(seedContent(startline + 6)) );
        if isempty(numIons) || abs(sum(floor(numIons))-sum(numIons)) > 0.001 
            seedProcessError(numOfPOSCAR, startline + 6);
            continue;
        end
        
        if ( isempty( strfind(char(seedContent(startline + 7)), 'D')) ...
                &&  isempty( strfind(char(seedContent(startline + 7)), 'd')) )
            seedProcessError(numOfPOSCAR, startline + 7);
            continue;
        end
        
        coord   = str2num( char(seedContent(startline + 8:startline+7+sum(numIons))) );
        if isempty(coord)
            seedProcessError(numOfPOSCAR, startline + 8);
            continue;
        end
        numOfPOSCAR = numOfPOSCAR +1;
        POSCARs(numOfPOSCAR).scale    = scale;
        POSCARs(numOfPOSCAR).atomType = atomType;
        POSCARs(numOfPOSCAR).numIons  = numIons;
        POSCARs(numOfPOSCAR).lattice  = lattice;
        POSCARs(numOfPOSCAR).coord    = coord;
        POSCARs(numOfPOSCAR).numMols  = numMols;
        
        %POSCARs(numOfPOSCAR)
        
        startline = startline+7+sum(numIons);
    catch
        warningStr = (['Seeds : Meet a problem when reading Seeds-' num2str(numOfPOSCAR),...
    ' starting at line ' num2str(startline) ', stop reading the POSCAR file']);
        USPEXmessage(0, warningStr, 0);
    end
end






function numMolsRead = numMolReadSeeds( string0 )


[nothing, numIons] = unix([' echo "' string0 '" |cut -d[  -f2 |cut -d] -f1  ']);
numMolsRead = str2num(numIons);



function seedProcessError( numSeed, startline )

warningStr = (['Seeds : Meet a problem when reading Seeds-' num2str(numSeed),...
    ' starting at line ' num2str(startline) ', trying to extract correct POSCARs']);
USPEXmessage(0, warningStr, 0);
%disp(warningStr);
