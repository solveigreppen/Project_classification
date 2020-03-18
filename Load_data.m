
%vowels= ['ae', 'ah', 'aw', 'eh', 'er', 'ei', 'ih', 'iy', 'oa', 'oo', 'uh', 'uw']; 


function table = load_data(filename)
if nargin < 1
    filename = 'vowdata_nohead.dat';
end

fileID = fopen(filename);
formatSpec = '%5s%4f%4f%5f%5f%5f%5f%5f%5f%5f%5f%5f%5f%5f%5f%5f%C%[^\n\r]';
C = textscan(fileID,formatSpec);

% Remove white space around all cell columns.
%dataArray{1} = strtrim(dataArray{1});

% Close the text file.
fclose(fileID);

end




