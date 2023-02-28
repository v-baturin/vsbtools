function [lat1,Surcandidate,surnumIons]=reduce_Surface(lattice,coordinates,numIons,Ind,varcomp)
%The new idea is implemented to improve the efficiency. 
%Sometimes, the surface atoms will get far away from the bulk, 
%and then some extreme horriable structures appear. 
%To avoid this, we firstly fix all the bulk atoms in the first step. 
%If the surface atoms are not bonded but go to the vacuum, just delete them.    
global ORG_STRUC 
global POP_STRUC

flag = 0;  %if rebuild POP
Step       = POP_STRUC.POPULATION(Ind).Step-1;
chanAList  = POP_STRUC.POPULATION(Ind).chanAList;
MaxStep    = length([ORG_STRUC.abinitioCode]);
thicknessS = ORG_STRUC.thicknessS;
        d1 = ORG_STRUC.bulk_lat(3,3);    %bulk length
        d2 = ORG_STRUC.thicknessS;       %surface
        d3 = ORG_STRUC.vacuumSize(Step); %vacuum
  atomType = ORG_STRUC.atomType;
  atomtype = [];
for i = 1:length(numIons)
    atomtype = [atomtype; atomType(i)*ones(numIons(i),1)];
end

Netcandidate=[];  %%%to store the add atoms in the network
NetatomType=[];   %%%
%%%%%%%%%%%%%%%%%%
if Step < MaxStep
    [isAllConnected, isConnectedToSubstrate ] = ...
    surface_connectivity_check( lattice, coordinates, atomtype, chanAList);

    if ~isAllConnected
       disp(['Structure ' num2str(Ind) ' not AllConected at Step ' num2str(Step)]);
       flag = 1;
    end
else
   isConnectedToSubstrate = ones(sum(numIons),1);
end

count = 0;
for i=1:sum(numIons)
   if chanAList(i) && isConnectedToSubstrate(i)
      count = count+1;
      Netcandidate(count,:) = coordinates(i,:);
      NetatomType(count)    = atomtype(i);
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Surcandidate=[];
SuratomType=[];
lat1= lattice;
lat1(3,3)= d2;
if ~isempty(Netcandidate)        
   Netcandidate(:,3)=Netcandidate(:,3)*lattice(3,3);
   if Step < MaxStep 
       item=0;
       for i=1:count
          if( Netcandidate(i,3) < d1+d2) && (Netcandidate(i,3)>d1-1)
             item=item+1;
             Surcandidate(item,:)=Netcandidate(i,:);
             SuratomType(item)=NetatomType(i);
          else
             flag=1; %rebuild
             disp(['delete the atoms going outside the surface']);
          end                    
       end
   else
       Surcandidate=Netcandidate;
       SuratomType=NetatomType;
   end

   for i=1:size(Surcandidate,1)
      Surcandidate(i,3) = (Surcandidate(i,3)-d1+0.1)/d2;
   end
end

surnumIons=[];
for i=1:length(ORG_STRUC.atomType)
     surnumIons(1,i)=sum( SuratomType==ORG_STRUC.atomType(i));    
end
%%%%%%%%%%%rebuild the POP_STRUC
if varcomp
   if (flag==1) && (Step < MaxStep)
      disp('rebuild the POP_STRUC');
      bulk_lattice=POP_STRUC.POPULATION(Ind).Bulk_LATTICE;
      bulk_pos=POP_STRUC.POPULATION(Ind).Bulk_COORDINATES;
      bulk_numIons=POP_STRUC.POPULATION(Ind).Bulk_numIons;
      [lat,coor,numIons,chanAList] = ...
      makeSurface(lat1,Surcandidate,surnumIons,bulk_lattice,bulk_pos, bulk_numIons, d3); 
      POP_STRUC.POPULATION(Ind).COORDINATES = coor;
      POP_STRUC.POPULATION(Ind).numIons = numIons;
      POP_STRUC.POPULATION(Ind).LATTICE = lat;
      POP_STRUC.POPULATION(Ind).chanAList=chanAList;
   end
end
