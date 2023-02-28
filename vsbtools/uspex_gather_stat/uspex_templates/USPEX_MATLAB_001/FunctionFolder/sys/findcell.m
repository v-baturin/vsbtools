function hk_family = find_cell(lat, Max_multi)
lat(3,:)=[];
lat(:,3)=[];
a = lat(1,:);
b = lat(2,:);
Max_iter = Max_multi;
phi = acos((dot(a,b))/norm(a)/norm(b));
hk = [];

for i = 1:Max_iter
   for j = 0:Max_iter
       for k = 0:Max_iter
           for l = 1:Max_iter
               a1 = i*a - j*b;
               b1 = k*a + l*b;
               multi = abs(det([a1;b1])/det(lat));
               if multi <= Max_multi
                  phi1 = acos((dot(a1,b1))/norm(a1)/norm(b1));
                  if abs(phi-phi1) < 0.001
                     hk = [hk; [i,-j,k,l]];
                  end
               end
           end
       end
   end
end

hk_family = hk(1,:);
for i = 2:size(hk,1)
    a1 = hk(i,1)*a + hk(i,2)*b;
    b1 = hk(i,3)*a + hk(i,4)*b;
    toAdd = 1;
    for j=1:i-1
        a2 = hk(j,1)*a + hk(j,2)*b;
        b2 = hk(j,3)*a + hk(j,4)*b;
        if (abs(norm(a2)/norm(a1)-1)<0.01) && (abs(norm(b2)/norm(b1)-1)<0.01)
           toAdd = 0;
           break;
        elseif (abs(norm(a2)/norm(b1)-1)<0.01) && (abs(norm(b2)/norm(a1)-1)<0.01)
           toAdd = 0;
           break;
        end
    end
    if toAdd
       hk_family = [hk_family; hk(i,:)];
    end
end


