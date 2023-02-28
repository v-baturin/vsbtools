function Clean_AbinitCode(code)

Files = {};
if code == 1
   movefile('CONTCAR', 'CONTCAR_old');
   movefile('OUTCAR',  'OUTCAR_old');
   movefile('OSZICAR', 'OSZICAR_old');
   movefile('POSCAR',  'POSCAR_old');
elseif code ==2
   Files = {'*-pbe', '*.ion*', 'siesta.*', 'INPUT_TMP.*', '*.psdump', '*.confpot'};
elseif code ==3
   Files = {'input', 'optimized.structure'};
   movefile('output', 'gulp-old.output');
elseif code ==6
   Files = {'fort.12', 'fort.13', 'dma_output', 'mol.dmain'};
elseif code ==7
   movefile('USPEX-1.cell', 'old_USPEX-1.cell');
   movefile('USPEX-pos-1.xyz', 'old_USPEX-pos-1.xyz');
   Files = {'USPEX*',  '*.uspex',  'cp2k.inp',  'lattice.stress_tensor'};
elseif code ==8
   movefile('output', 'qE_old.output');
elseif code ==9
   movefile('FHI_output', 'FHI_output_old');
   Files = {'relaxation_restart_file.FHIaims'};
elseif code == 12
   Files = {'output', '*.current_stage', '*.int*', '*.key', '*.make0', ...
  '*.make','*.seq', '*.tmp', '*.xyz*', '*.angles*', '*.pdb*', 'POSCAR*'};
elseif code == 14
  Files = {'BoltzTraP.def', 'BoltzTraP.output', 'CalcFold*.condtens*', ...
           'CalcFold*.en*', 'CalcFold*.hall*', 'CalcFold*.intrans', ...
           'CalcFold*.lfengre', 'CalcFold*.outputtrans', 'CalcFold*.ph*', ...
           'CalcFold*.sig*', 'CalcFold*.struct', 'CalcFold*.trace*', ...
           'CalcFold*.transdos', 'CalcFold*.v2dos' };
elseif code == 15
   Files = {'detailed.out'};
elseif code ==17
	Files = {'BLANDAD', 'CHG','DISP','DOS*','EIGENVAL', 'F*', 'HARMONIC', ...
	         'IBZKPT', 'KRAFTER*','MEMORYFILE*',' OUTCAR*', ' POSCAR.phon', 'POSCAR_REF','POSCAR_START*', ...
		  'POSCA*','PROJECTIONS.dat','QPOINTS','SCPH_IS_DONE',' SPOSCAR', 'SYMOP', 'TRANSROT', ...
		   'vasprun.xml', 'WAVECAR','CONVERGENCE' };
end

DeleteFiles(Files);
