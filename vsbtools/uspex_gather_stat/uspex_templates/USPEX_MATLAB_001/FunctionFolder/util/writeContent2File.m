function writeContent2File(file, content, option)


f = fopen(file, option);

for i=1:size(content, 1)
    fprintf(f, '%s\n', char(content(i)));
end

fclose(f);
