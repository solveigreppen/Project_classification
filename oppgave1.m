[files, dur, F0s, F1s, F2s, F3s, F4s, F120, F220, F320, F150, F250, F350, F180, F280, F380] = textscan('vowdata_nohead.dat', '%s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f');

vowels= char('ae', 'ah', 'aw', 'eh', 'er', 'ei', 'ih', 'iy', 'oa', 'oo', 'uh', 'uw'); 
talker_group = char('m', 'w', 'b', 'g'); 
filenames = char(filenames); %converts cell array to character matrix.
[nfiles, nchar] = size(filenames); 

for ifile=1:nfiles,
    vowel_code(ifile) strmatch(filenames(ifile, 4:5), vowel); 
    
    





% % character 1:     m=man, w=woman, b=boy, g=girl 
% % characters 2-3:  talker number 
% % characters 4-5:  vowel (ae=”had”, ah=”hod”, aw=”hawed”, eh=”head”,
% %    er=”heard”, ei=”haid”, ih=”hid”, iy=”heed”, oa=/o/ as in “boat”,
% %    oo=”hood”, uh=”hud”, uw=”who’d”)
% 
% [files, dur, F0s, F1s, F2s, F3s, F4s, F120, F220, F320, F150, F250, F350, F180, F280, F380] = 