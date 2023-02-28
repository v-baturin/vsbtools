function sQE = structureQuasiEntropy(numIons, atom_fing, weight)

sQE = 0;
[trash, r, c] = size(atom_fing);
for i = 1 : length(numIons)
   if numIons(i) > 1
      tmp = 0;
      k = 0;
      start_n = sum(numIons(1:i-1));
      for j = start_n+1 : start_n+numIons(i)-1
          tmp_fing1 = reshape(atom_fing(j,:,:),r,c);
          for j1 = j+1 : start_n+numIons(i)
            k = k + 1;
            tmp_fing2 = reshape(atom_fing(j1,:,:),r,c);
            dist = cosineDistance(tmp_fing1, tmp_fing2, weight);
            tmp = tmp + (1 - dist)*log(1 - dist);
          end
      end
      sQE = sQE + numIons(i)/sum(numIons)*tmp/k;
   end
end
sQE = -sQE;


