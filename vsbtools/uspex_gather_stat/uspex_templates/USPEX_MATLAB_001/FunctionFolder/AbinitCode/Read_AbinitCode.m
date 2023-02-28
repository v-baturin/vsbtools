function target = Read_AbinitCode(code, flag, ID, Step)
% 0: if it is complete
%-1: if 1st SCF is done or if SCF is converged
% 1: energy/enthalpy
% 2: pressue tensor
% 3: dielectric constant
% 4: gap
% 5: magmoment
%-5: atomic magnetization 
% 6: elestic properties

if code == 1
   target = Read_VASP(flag, ID, Step);  %current (-1, 0, 1, 2, 3, 4, 5, -5, 6)
elseif code ==2
   target = Read_SIESTA(flag); %current (-1, 0, 1)
elseif code ==3
   target = Read_GULP(flag, ID);  % current (-1, 0, 1, 2, 3, 6)
elseif code ==4
   target = Read_LAMMPS(flag); % current (-1, 0, 1)
elseif code ==5
   target = Read_NeurNetw(flag); % current (-1, 0, 1)
elseif code ==6
   target = Read_DMACRYS(flag); % current ( -1, 0, 1)
elseif code ==7
   target = Read_CP2K(flag); % current (-1, 0, 1)
elseif code ==8
   target = Read_QE(flag, ID); % current (-1, 0, 1)
elseif code ==9
   target = Read_FHIaims(flag, ID); % current (-1, 0, 1)
elseif code ==10
   target = Read_ATK(flag); % current (-1, 0, 1)
elseif code ==11
   target = Read_CASTEP(flag); % current (0, 1)
elseif code == 12
   target = Read_Tinker(flag); % current (0, 1)
elseif code == 13
   target = Read_MOPAC(flag, ID); % current (-1, 0, 1)
elseif code == 14 %% thermoelectric properties
    if flag == 1
        target = Read_VASP(flag, ID);   % Read enthalpy
    else
        % flag == 2 power factor, deprecated
	% flag == 3 ZT, thermoelectric property
        target = Read_BoltzTraP(flag, ID); 
    end
elseif code == 15 %% DFTB
   target = Read_DFTB(flag, ID); % current (-1, 0, 1)
elseif code == 16
    target = Read_Gaussian(flag, ID);
elseif code ==17
	if flag ==1
		target = Read_VASP(flag, ID); %include Fel(V,T) term by using Fermi smearing function
	else
	       target = Read_SCPH(flag, ID); %get Fphon from SCAILD calculation
  end
end
