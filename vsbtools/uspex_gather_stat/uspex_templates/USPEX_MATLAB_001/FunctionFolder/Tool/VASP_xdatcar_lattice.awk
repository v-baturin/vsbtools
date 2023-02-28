BEGIN { state="header" }
state == "header"  { getline; state="lattice"; }
state == "lattice" { for(i=1; i<=3; i++) {
                       getline; print $1 "," $2 "," $3;
                     }
                     state="atoms"
                   }
state == "atoms"   { getline; getline;
                     natoms=0; for(i=1; i<=$NF; i++) { natoms+=$i; }
                     getline;
                     state="info" }
state == "info"    { if($0 ~ /configuration/) {
                       state = "coords";
                     } else {
                       state = "lattice";
                     }
                   }
state == "coords"  { for(i=1; i<=natoms; i++) { getline; } state="info" }
