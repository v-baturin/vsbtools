function TPSAbinitWriting(code, Ind)



if     code ==4
    Write_LAMMPS_TPS(Ind)
elseif code ==7
    Write_CP2K_TPS(Ind)
end



%
%
%
function     Write_CP2K_TPS(Ind)
Ind;
