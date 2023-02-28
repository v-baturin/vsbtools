function DeleteFiles(Files)

%Octave (up to 3.8) does not support Delete Mutiple Files at once

for i = 1:length(Files)
    delete(Files{i});
end
