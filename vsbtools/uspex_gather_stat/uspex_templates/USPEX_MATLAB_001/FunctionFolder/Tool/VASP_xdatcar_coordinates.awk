BEGIN { state="header" }
state == "header"  { getline; getline; state="lattice"; }
state == "lattice" { for(i=1; i<=3; i++){ getline } state="atoms" }
state == "atoms"   { getline;
                     natoms=0; for(i=1; i<=$NF; i++) { natoms+=$i };
                     getline;
                     state="info" }
state == "info"    { if($0 ~ /configuration/) {
                       state = "coords";
                     } else {
                       getline;
                       state = "lattice";
                     }
                   }
state == "coords"  { for(i=1; i<=natoms; i++) {
                     getline; print $1 "," $2 "," $3
                   }
                   state="info"}
