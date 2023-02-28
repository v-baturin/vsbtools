function [candidate, Lattice, numSites, Operation, errorMsg] = ...
                                  symope_3D_MOL(nsym, numIons, lat, minD)

global ORG_STRUC

fixLat = ORG_STRUC.constLattice;        %-[QZ]We must need to know from ORG_STRUC
spgBINDIR=[ORG_STRUC.USPEXPath '/FunctionFolder/spacegroup']; % path of execution
fixRndSeed = ORG_STRUC.fixRndSeed;          % fix rand seed, to reproduce results
sGroup = spaceGroups(nsym);                       % space group's standard symbol

%----------Flowchart-------------
if fixLat
  [P, PB] = GetPermutation(nsym);                           %1, permutation axis
else
  P  = [1 2 3];
  PB = [1 2 3];
end
Init_lat   = Get_Init_Lattice(lat, P);                    %2, Get initial lattice
[Init_numIons, Init_lat,fail] = GetPrimitiveCell(sGroup, Init_lat, numIons,fixLat);
                                                          %4, Primitive cell
if ~fail
                                                 %5, Prepare input for Stokes code
    Write_Stokes_input(Init_numIons, minD, nsym, Init_lat, fixRndSeed, []); 
    unixCmd([spgBINDIR '/random_cell_mol < rc.in > rc.out']);     %6, Execute code
    [Operation, numSites, candidate, Lattice, errorMsg] = ...
                       Read_Stokes_output_MOL('rc.out');   %7, Read Stokes outputs
    if errorMsg ==0                                  %-8, Adjust lattice if needed
       Lattice(4:6) = Lattice(4:6)*pi/180;
       Lattice = latConverter(Lattice);
       %disp('primitive:')
       %disp(num2str(candidate))
       if (fixLat) && (sGroup(1) ~= 'P') 
          [Lattice, candidate, numIons0, Operation, errorMsg] = find_symmetry...
                                     (Lattice, candidate, Init_numIons, 1, 0.01 );
          if errorMsg > 0
              disp(['error in finding symmetry: ' num2str(nsym)]);
          else
              duplicate = sum(numIons0)/sum(Init_numIons);   %duplication
              numSites = numSites*duplicate;
              numIons =  Init_numIons*duplicate;
          end
       end
        
       %disp('Conventional:')
       %disp(num2str(candidate))
       if size(candidate, 1) ~= sum(numIons)
          errorMsg = 1;
       else
          [Lattice, candidate] = Get_Final_Struc(Lattice, candidate, PB);
          Operation = Permutation_symmetry(Operation, PB);

          %disp('TTTTTTTTTTTTTTTTTTTTTTTT')
          %Lattice
          %candidate
          %P
       end
    end
end


function Operation = Permutation_symmetry(Operation, PB)
%Here we only need to consider Pmm2, so just swap things like
%( x  y  z)    x     ( x  y  z)
%(-x  y  z)  ---->   ( x  y -z) 
%( x -y  z)  <----   ( x -y  z) 
%(-x  y  z)    z     ( x  y -z)

%must be diagonal matrix, we only need to change the sign.
N_sym = size(Operation,1)/4;
for i = 1:N_sym
    Opt = Operation( (i-1)*4+1 : i*4, : );
    Opt1 = Opt;
    for j = 1:3
        Opt(j,j) = Opt1(PB(j),PB(j));
        Opt(4,j) = Opt1(4,PB(j));
    end
    Operation( (i-1)*4+1 : i*4, : ) = Opt;
end
