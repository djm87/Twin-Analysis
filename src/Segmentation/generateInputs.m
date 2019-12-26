function [outputArg1,outputArg2] = generateInputs(flist)
%generateInputs adds files to directory if they don't exist assuming flist is 
%a cell array of strings
for i=1:length(flist)
    if ~isfile(flist{i})
        fid = fopen(flist{i},'w');
        fclose(fid);
    end
end
end

