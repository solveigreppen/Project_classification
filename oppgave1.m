[files, dur, F0s, F1s, F2s, F3s, F4s, F120, F220, F320, F150, F250, F350, F180, F280, F380] = textread('vowdata_nohead.dat', '%s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f');

vowels= ['ae'; 'ah'; 'aw'; 'eh'; 'er'; 'ei'; 'ih'; 'iy'; 'oa'; 'oo'; 'uh'; 'uw'];
talker_group = char('m', 'w', 'b', 'g'); 
filenames = char(files); %converts cell array to character matrix.
[nfiles, nchar] = size(filenames); %finner ut hvor bred og lang matrisen er(?) n files er nedover
len_vowels = length(vowels); 
% for ifile=1:nfiles
for vowel=1:len_vowels
    a = strcmp(filenames(5, 4:5), vowels(vowel, 1:2)); %filenames(ifile, 4:5) gives us the row, and our vowel, that is in place 4 and 5
    if a=1 
        
disp(a);
disp(vowels(1, 1:2));
disp(filenames(5, 4:5));
s1= 'kajak'; 
s2= 'kaja'; 
tf=strncmp(s1,s2,4);
disp(tf); 

    
    





% % character 1:     m=man, w=woman, b=boy, g=girl 
% % characters 2-3:  talker number 
% % characters 4-5:  vowel (ae=”had”, ah=”hod”, aw=”hawed”, eh=”head”,
% %    er=”heard”, ei=”haid”, ih=”hid”, iy=”heed”, oa=/o/ as in “boat”,
% %    oo=”hood”, uh=”hud”, uw=”who’d”)
% 
% [files, dur, F0s, F1s, F2s, F3s, F4s, F120, F220, F320, F150, F250, F350, F180, F280, F380] = 