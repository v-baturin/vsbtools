function [candidate, Lattice, errorS] = symope_crystal(nsym, numIons, lat,...
                                                                minD, sym_coef)

global ORG_STRUC

fixLat = ORG_STRUC.constLattice;    %-[QZ]We must need to know from ORG_STRUC

spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup']; % path of execution
fixRndSeed = ORG_STRUC.fixRndSeed;          % fix rand seed, to reproduce results
sGroup = spaceGroups(nsym);                 % space group's standard symbol

% ---------Initilizating outputs------------
candidate = zeros(sum(numIons),3);                       % 1, init_coordinates
Lattice = [1 0 0; 0 1 0; 0 0 1];                         % 2, init_lattice
errorS = 0;                                              % 3, error check

%----------Flowchart-------------
[permutation, permutationBack] = GetPermutation(nsym);    %1, permutation axis
Init_lat   = Get_Init_Lattice(lat, permutation);          %2, Get initial lattice
%Init_lat = fix_latticeStokes_before(nsym, Init_lat);     %3, fix non conventional                                            
                                           %choice of primitive monoclinic lattice
[Init_numIons, Init_lat,fail] = GetPrimitiveCell(sGroup, Init_lat, numIons,fixLat);     
                                                                %4, Primitive cell

if ~fail
                                                 %5, Prepare input for Stokes code
    Write_Stokes_input(Init_numIons, minD, nsym, Init_lat, fixRndSeed, sym_coef); 
    unixCmd([spgBINDIR '/random_cell < rc.in > rc.out']);  %6, Execute Stokes code
    [coordinate_S, lattice_S, failed]=Read_Stokes_output('rc.out');               
                                                           %7, Read Stokes outputs
    
    if ~failed                                      %-8, Adjust lattice if needed
      [lattice, coordinate]=fix_latticeStokes_after(nsym, lattice_S, coordinate_S);
      lattice(4:6) = lattice(4:6)*pi/180;
      if fixLat == 1
          Lattice_Matrix = latConverter(lattice);
      else
          Lattice_Matrix = latConverter(lattice);
          abs_cand = coordinate*Lattice_Matrix;
          [abs_cand,Lattice_Matrix] = optLattice(abs_cand, Lattice_Matrix);     
                                                    %--------8.1  optimize lattice
          coordinate = abs_cand/Lattice_Matrix;
      end
      Lattice_Matrix(find(abs(Lattice_Matrix)<0.000001))=0;
      coordinate(find(abs(coordinate  )<0.000001))=0;
      coordinate(find(abs(coordinate-1)<0.000001))=0;
        
        if latticeCheck(Lattice_Matrix)             %-9, GetConventional cell 
            Lat_type = sGroup(1);
            if (Lat_type ~= 'P') && (fixLat == 1)
                if (nsym > 37) && (nsym < 42)        %-9.1 special  groups
                    Lat_type = 'B';
                end
                [Lattice_Matrix, candidate] = unitCellFromPrimitive...
                                (Lat_type, Lattice_Matrix, coordinate);
            else
                candidate = coordinate;
            end                                     %-10  Permuation back
            [Lattice, candidate] = Get_Final_Struc...
                            (Lattice_Matrix, candidate, permutationBack);
        else
            errorS = 1;
        end
        % ---------Fix for ticket #152 - Problem with splitBigCell 
        if size(candidate, 1) ~= sum(numIons)
            errorS = 1;
        end
    else
        errorS = 1;
    end
end
