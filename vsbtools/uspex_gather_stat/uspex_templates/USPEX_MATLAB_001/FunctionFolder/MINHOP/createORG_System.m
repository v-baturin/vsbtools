function createORG_System(inputFile)

global ORG_STRUC

getPy=[ORG_STRUC.USPEXPath,'/FunctionFolder/getInput.py'];

%%%%%CalcType
%3 digits (dimension: 0-3; molecule: 0/1; varcomp: 0/1)
calculationType = python_uspex(getPy, ['-f ' inputFile ' -b calculationType -c 1']);
if ~isempty(calculationType)
ORG_STRUC.dimension = str2num(calculationType(1));
ORG_STRUC.molecule  = str2num(calculationType(2));
ORG_STRUC.varcomp   = str2num(calculationType(3));
end

%pickUpYN = python_uspex(getPy, ['-f ' inputFile ' -b pickUpYN -c 1']);
%if isempty(pickUpYN)
%  pickUpYN = '0'; % default
%end
%ORG_STRUC.pickUpYN = str2num(pickUpYN);

ORG_STRUC.pickedUP=0;
pickUpGen = python_uspex(getPy, ['-f ' inputFile ' -b pickUpGen -c 1']);
if isempty(pickUpGen)
  pickUpGen = '1'; % default
else
  if str2num(pickUpGen) > 0
     ORG_STRUC.pickUpYN = 1;
  end
end

ORG_STRUC.pickUpGen = str2num(pickUpGen);
if ORG_STRUC.pickUpYN == 0
 ORG_STRUC.pickUpGen = 1; % needed for BestEnthalpies plot
end

pickUpFolder = python_uspex(getPy, ['-f ' inputFile ' -b pickUpFolder -c 1']);
if isempty(pickUpFolder)
  pickUpFolder = '1'; % default
end
ORG_STRUC.pickUpFolder = str2num(pickUpFolder);
%%%%%%%read initial structure%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fid,message] = fopen('POSCAR_1');
tmp = fgetl(fid) ;% system description
scale_factor = fgetl(fid) ;% 1 = numbers in angstrems
lat = fscanf(fid,'%g',[3,3]); % lattice vectors
lat = lat';
tmp = fgetl(fid);
atomType = fgetl(fid);
numIons = fgetl(fid);
if strcmp(numIons, 'Direct') % VASP 4
  ORG_STRUC.numIons = str2num(atomType);
 % types of ions; can be number, short name of full name. USPEX will use the numbers only
 atomType = python_uspex(getPy, ['-f ' inputFile ' -b atomType -e EndAtomType']);
 atomType(end) = [];
else    % VASP5
  ORG_STRUC.numIons = str2num(numIons);
  tmp = fgetl(fid);
end

%ORG_STRUC.atomType = str2num(atomType);
ORG_STRUC.atomType = zeros(1,size(ORG_STRUC.numIons,2));
c1 = findstr(atomType, ' ');
c = sort(str2num(['0 ' num2str(c1)]));
c(end+1) = length(atomType) + 1;
ind1 = 1;
for i = 2 : length(c)
 if c(i-1)+1 > c(i)-1
   continue
 end
 tmp = atomType(c(i-1)+1 : c(i)-1);
 if ~isempty(str2num(tmp)) % number
   ORG_STRUC.atomType(ind1) = str2num(tmp);
 else
   for j = 1 : 105
    if strcmp(lower(tmp), lower(elementFullName(j))) || strcmp(lower(tmp), lower(megaDoof(j)))
      ORG_STRUC.atomType(ind1) = j;
      break;
    end
   end
 end
 ind1 = ind1 + 1;
end

sss = fscanf(fid,'%g',[3,sum(ORG_STRUC.numIons)]);
ss = sss';
coordinates = ss(:,1:3);
     fclose(fid);
%%%%%%%lattice conversion%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base(1,:) = lat(1,:)/norm(lat(1,:));
base(3,:) = cross(lat(1,:),lat(2,:));
base(3,:) = base(3,:) / norm(base(3,:));
base(2,:) = cross(base(3,:),base(1,:));
ORG_STRUC.lattice = lat/base;
ORG_STRUC.coordinates = coordinates - floor(coordinates);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% BEGIN LATTICE manipulations %%%%%%%%%%%%
 % varcomp does not support fixed lattice
 lattice_INPUT = python_uspex(getPy, ['-f ' inputFile ' -b Latticevalues -e Endvalues'], 1);
 
 if size(lattice_INPUT, 1) == 1
     if (size(lattice_INPUT, 2) == 6)     % fixed lattice in a form {a b c alpha beta gamma}?
         lattice_INPUT(4:6) = lattice_INPUT(4:6)*pi/180;         % Matlab works with radians, input - in degrees
         ORG_STRUC.constLattice = 1;
         ORG_STRUC.lattice = latConverter(lattice_INPUT);   % 3x3 form
         ORG_STRUC.latVolume = det(ORG_STRUC.lattice);
     else  % varcomp
         ORG_STRUC.latVolume = lattice_INPUT;
         ORG_STRUC.constLattice = 0;
     end
 elseif (sum(size(lattice_INPUT) == [3 3]) == 2)  %3*3 matrix
     ORG_STRUC.latVolume = det(lattice_INPUT);
     ORG_STRUC.constLattice = 1;
     ORG_STRUC.lattice = lattice_INPUT;
 else %empty
 %    if ORG_STRUC.dimension~=2 & ORG_STRUC.dimension~=-3 & ORG_STRUC.dimension~=-4
 %       estimateInputVolume();
        ORG_STRUC.constLattice = 0;
 %    end
 end
 
 if ORG_STRUC.dimension==0 || ORG_STRUC.dimension==2 || ORG_STRUC.dimension==-3
    ORG_STRUC.constLattice = 1;
 end
 %%%%%%% END LATTICE manipulations %%%%%%%%%%%%


ORG_STRUC.goodBonds = zeros(length(ORG_STRUC.atomType));
  for i = 1 : length(ORG_STRUC.atomType)
   for j = i : length(ORG_STRUC.atomType)
     ORG_STRUC.goodBonds(i,j) = 0.15;
     ORG_STRUC.goodBonds(j,i) = 0.15;
   end
  end
% valences for each type of atoms
valences = python_uspex(getPy, ['-f ' inputFile ' -b valences -e endValences']);
ORG_STRUC.valences = str2num(valences);
if isempty(valences) % default
  ORG_STRUC.valences = zeros(1,length(ORG_STRUC.atomType));
  for i = 1 : length(ORG_STRUC.atomType)
      ORG_STRUC.valences(i) = str2num(valence(ORG_STRUC.atomType(i)));
  end
else
  ORG_STRUC.valences = str2num(valences);
end

NvalElectrons = python_uspex(getPy, ['-f ' inputFile ' -b valenceElectr -e endValenceElectr']);
if isempty(NvalElectrons) % default
  ORG_STRUC.NvalElectrons = zeros(1,length(ORG_STRUC.atomType));
  for i = 1 : length(ORG_STRUC.atomType)
    ORG_STRUC.NvalElectrons(i) = str2num(valenceElectronsNumber(ORG_STRUC.atomType(i)));
    if ORG_STRUC.NvalElectrons(i) == 0  % transition metal (d-, f-)
      ORG_STRUC.NvalElectrons(i) = ORG_STRUC.valences(i);
    end
  end
else
  ORG_STRUC.NvalElectrons = str2num(NvalElectrons);
end

hardCore = python_uspex(getPy, ['-f ' inputFile ' -b IonDistances -e EndDistances']);
if isempty(hardCore)
  ORG_STRUC.minDistMatrice = zeros(1,length(ORG_STRUC.atomType));
  for i = 1 : length(ORG_STRUC.atomType)
    s = covalentRadius(ORG_STRUC.atomType(i));
    ORG_STRUC.hardCore(i) = str2num(s)/2;
  end
else
  ORG_STRUC.hardCore = str2num(hardCore);
end

%%%% check whether hard core radii or distance matrix were given and create matrix if necessary %%%%%%%
ORG_STRUC.minDistMatrice = zeros(length(ORG_STRUC.atomType), length(ORG_STRUC.atomType));
if size(ORG_STRUC.hardCore,1) == size(ORG_STRUC.hardCore,2)
  for i = 1 : length(ORG_STRUC.atomType)
   for j = i : length(ORG_STRUC.atomType)
     ORG_STRUC.minDistMatrice(i,j) = ORG_STRUC.hardCore(i,j);
     ORG_STRUC.minDistMatrice(j,i) = ORG_STRUC.hardCore(i,j);
   end
  end
else
  ORG_STRUC.minDistMatrice = zeros(length(ORG_STRUC.hardCore));
  for hardInd_1 = 1 : length(ORG_STRUC.hardCore)
    for hardInd_2 = hardInd_1 : length(ORG_STRUC.hardCore)
       ORG_STRUC.minDistMatrice(hardInd_1,hardInd_2) = ORG_STRUC.hardCore(hardInd_1) + ORG_STRUC.hardCore(hardInd_2);
       ORG_STRUC.minDistMatrice(hardInd_2,hardInd_1) = ORG_STRUC.hardCore(hardInd_1) + ORG_STRUC.hardCore(hardInd_2);
    end
  end
end

ORG_STRUC.bestBasicStructure = ORG_STRUC.coordinates;
