gen=`tail -1 results1/BESTIndividuals | sed 's/^ *//;s/ *$//' | cut -d ' ' -f 1`
enth=`tail -1 results1/BESTIndividuals | grep -o -P '(?<=\]).*(?=\[)' | sed 's/^ *//;s/ *$//' | cut -d ' ' -f 1`
echo $gen $enth
