function WriteTEproperties(resFolder)
% Write the thermoelectric properties to files.
% Current implementation is based on file renaming.
% TE data are too large to be stored in USPEX.mat.

global POP_STRUC;

% variables
gen = POP_STRUC.generation;
order = [POP_STRUC.POPULATION(:).Number];
suffixes = {'ZT.txt', 'ZT_xx.txt', 'ZT_yy.txt', 'ZT_zz.txt'};

% Folders and files
folder = [resFolder '/TEproperties'];
summary_file = [folder '/summary.txt'];

[order_sorted, index] = sort(order);
for jj = 1:length(order)
  Ind_No = index(jj);
  number = POP_STRUC.POPULATION(Ind_No).Number;
  file_prefix = ['Gen' num2str(gen) '-Ind' num2str(Ind_No) '-Step*'];

  summary_lines = dir([folder '/' file_prefix 'summary.txt']);
  if length(summary_lines) > 0
    fn = [folder '/' summary_lines(1).name];
    awk_cmd = ['awk ''{printf "%5d %s\n", ' num2str(number) ', $0}'''];
    unixCmd([awk_cmd ' < ' fn ' >> ' summary_file]);
    unixCmd(['rm ' fn]);
  end

  for kk = 1:length(suffixes)
    BoltzTraP_files = dir([folder '/' file_prefix suffixes{kk}]);
    if length(BoltzTraP_files) > 0
      src_fn = [folder '/' BoltzTraP_files(1).name];
      dst_fn = [folder '/' num2str(number, '%05d') '-' suffixes{kk}];
      unixCmd(['mv ' src_fn ' ' dst_fn]);
    end
  end
end
