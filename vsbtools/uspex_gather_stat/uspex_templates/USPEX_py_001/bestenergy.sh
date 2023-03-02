gen=`grep 'Generation' results1/BESTIndividuals | tail -1 | sed 's/  */ /g' | cut -d' ' -f 2`
enth=`tail -2 results1/BESTIndividuals | head -1 | sed 's/  */ /g' | cut -d' ' -f 9`
echo $gen $enth
