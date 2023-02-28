function Molecules = GetOrientation(candidate, lat, numSites, Operation, MtypeLIST, nsym)
global ORG_STRUC
Molecules=struct('MOLCOORS',{},'ZMATRIX',{},'ID',{},'MOLCENTER',{},'Operation',{});

Molecules=[];
STDMOL = ORG_STRUC.STDMOL;
%Molecules(1).MOLCOORS=STDMOL(MtypeLIST(1)).molecule;
num = size(Operation,1)/4;
goodcandidate = [];
good = 0;
lat_6 =  latConverter(lat);
lat_60 = latConverter([1 0 0; 0 1 0; 0 0 1]);

count = 0;
for i = 1:length(numSites)
    num_opt      = STDMOL(MtypeLIST(count+1)).num_optFlags;
    format       = STDMOL(MtypeLIST(count+1)).format;
    coords       = STDMOL(MtypeLIST(count+1)).molecule;
    coords = bsxfun(@minus, coords, mean(coords)); % Matrix-Vector;
    coords = Rotate_rigid_body([0 0 0],[1 0 0], coords, (rand-0.5)*2*pi);
    coords = Rotate_rigid_body([0 0 0],[0 1 0], coords, (rand-0.5)*2*pi);
    coords = Rotate_rigid_body([0 0 0],[0 0 1], coords, (rand-0.5)*2*pi);
    [Z, coords] = RotInertia(coords, format, num_opt, MtypeLIST(count+1));
    coords = bsxfun(@plus, coords, Frac2Cart(candidate(count+1,:), lat_6)); 

    for j=1:numSites(i)
        Opt = Operation( (j-1)*4+1 : j*4, : );
        R = Operate_molecule(coords, Opt, lat_6);
        Molecules(count+j).MOLCOORS= R;
        Molecules(count+j).ZMATRIX   = NEW_coord2Zmatrix(R,format);
        Molecules(count+j).Operation = Opt;
    end
    count = count + numSites(i);
end

