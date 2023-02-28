function Write_DFTB_gen(atomType, numIons, lattice, coor)

fp = fopen('DFTB.gen', 'w');
fprintf(fp, '%4d F\n', sum(numIons));
coor = coor-floor(coor);

if sum(atomType)>0
   for i=1:length(numIons)
       if numIons(i) > 0
          fprintf(fp,'%4s', megaDoof(ceil(atomType(i))));
       end
   end
   fprintf(fp, '\n');
end
count = 0;
for i = 1 : length(numIons)
    for j = 1:numIons(i)
        count = count + 1;
        fprintf(fp, '%4d %4d %12.6f %12.6f %12.6f\n', count, i, coor(count,:));
    end
end
fprintf(fp, '0.0000 0.0000 0.0000\n');
for latticeLoop = 1 : 3
   fprintf(fp, '%12.6f %12.6f %12.6f\n', lattice(latticeLoop,:));
end

fclose(fp);
