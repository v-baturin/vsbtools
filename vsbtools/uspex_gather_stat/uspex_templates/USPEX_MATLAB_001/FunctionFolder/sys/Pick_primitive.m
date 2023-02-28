function [newlattice,newcoord,numIons] = Pick_primitive(lattice, coord, numIons, lat_base, atomType)

newlattice=lattice;
newcoord=[];
newtypes=[];
type = [];
for i = 1:length(numIons)
   type = [type; atomType(i)*ones(numIons(i),1)];
end
num1 = floor(lattice(1,1)/lat_base(1,1)+0.01);
num2 = floor(lattice(2,2)/lat_base(2,2)+0.01);
%%%%%%%%%%%%%%how many atoms in each region
pick = [1,1];
max_atom = 0;
for i=1:num1
    for j=1:num2
        ioncount=0;
        X_max =     i/num1 - 0.0001;
        Y_max =     j/num2 - 0.0001;
        X_min = (i-1)/num1 - 0.0001;
        Y_min = (j-1)/num2 - 0.0001;
        %[X_min, X_max, Y_min, Y_max]
        for k=1:sum(numIons)
            if     ( X_min < coord(k,1) ) && ( coord(k,1) < X_max ) ...
                && (  Y_min< coord(k,2) ) && ( coord(k,2) < Y_max )
                ioncount=ioncount+1;
            end
        end
        if ioncount>max_atom
           max_atom = ioncount;
           pick = [i,j];
        end      
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
newlattice = lat_base;
X_max =     pick(1)/num1 - 0.0001;
Y_max =     pick(2)/num2 - 0.0001;
X_min = (pick(1)-1)/num1 - 0.0001;
Y_min = (pick(2)-1)/num2 - 0.0001;
ioncount = 0;

for i=1:sum(numIons)
    if     ( X_min < coord(i,1) ) && ( coord(i,1) < X_max ) ...
       && (  Y_min < coord(i,2) ) && ( coord(i,2) < Y_max )
        ioncount=ioncount+1;
        newcoord(ioncount,:)= coord(i,:)*lattice/newlattice;
        newtypes(ioncount)  = type(i);
    end
end

for i=1:size(atomType,2)
    numIons(i)=sum(size(find(newtypes==atomType(i)),2));
end
newcoord  = newcoord-floor(newcoord);
