
vowels= ['ae', 'ah', 'aw', 'eh', 'er', 'ei', 'ih', 'iy', 'oa', 'oo', 'uh', 'uw']; 


function table = load_vowdata(filename)
if nargin < 1
    filename = 'vowdata_nohead.dat';
end

