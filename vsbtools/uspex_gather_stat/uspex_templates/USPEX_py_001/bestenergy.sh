gen=`grep 'Generation' results1/BESTIndividuals | tail -1 | sed 's/  */ /g' | cut -d' ' -f 2`
enth=`tail -2 BESTIndividuals | head -1 | sed 's/  */ /g' | cut -d '|' -f 5`
echo $gen $enth
