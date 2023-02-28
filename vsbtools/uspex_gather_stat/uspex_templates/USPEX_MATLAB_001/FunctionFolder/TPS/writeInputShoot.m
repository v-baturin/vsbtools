function writeInputShoot()


global ORG_STRUC
global TPS_STRUC


% input_shoot file
numAtomType = size(ORG_STRUC.numIons,1);

amplitudeA2B = TPS_STRUC.amplitudeA2B;
amplitudeB2A = TPS_STRUC.amplitudeB2A;
magnitudeA2B = TPS_STRUC.magnitudeA2B;
magnitudeB2A = TPS_STRUC.magnitudeB2A;

fp = fopen('input_shoot','w');

fprintf(fp,'%d\n', numAtomType);

for i = 1:numAtomType
    whichAtom = find(ORG_STRUC.numIons(i,:)>0);
    if length( whichAtom ) == 1
        fprintf(fp,'%5d  %8.4f   %12s\n',  ORG_STRUC.numIons(i,whichAtom), ...
            ORG_STRUC.mass(i), ORG_STRUC.speciesSymbol{i} );
    else
        fprintf(fp,'%5d  %8.4f   %12s\n',  1, ORG_STRUC.mass(i), ORG_STRUC.speciesSymbol{i} );
    end
end

fprintf(fp, '%10.8f %10.8f \n', amplitudeA2B, magnitudeA2B);
fprintf(fp, '%10.8f %10.8f \n', amplitudeB2A, magnitudeB2A);

fprintf(fp, '%d\n', TPS_STRUC.iteration);

fclose(fp);


% lastrun
if TPS_STRUC.lastrun.success == 0
    success = -1;
else
    success =  1;
end
fp = fopen('lastrun','w');
fprintf(fp, '%d %s', success, TPS_STRUC.lastrun.direction(end));
fclose(fp);
