clear;
Ntrain = 70;
Ntest = 69;
[files,dur,F0s,F1s,F2s,F3s,F4s,F120,F220,F320,F150,F250,F350,F180,F280,F380] = textread('vowdata_nohead.dat','%s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f'); 
nfiles = 1668;
% character 1:     m=man, w=woman, b=boy, g=girl 
% characters 2-3:  talker number 
% characters 4-5:  vowel (ae=”had”, ah=”hod”, aw=”hawed”, eh=”head”,  
%    er=”heard”, ei=”haid”, ih=”hid”, iy=”heed”, oa=/o/ as in “boat”,  
%    oo=”hood”, uh=”hud”, uw=”who’d”) 


vowel = str2mat('ae','ah','aw','eh','er','ei','ih','iy','oa','oo','uh','uw'); 
talker_group = char('m','w','b','g');  

vowel_code = zeros(1,nfiles);
talker_number = zeros(1,nfiles);

files=char(files); % convert cell array to character matrix [nfiles,nchar]=size(filenames); 
for i=1:nfiles  
    vowel_code(i) = strmatch(files(i,4:5),vowel);  
    talker_group_code(i) = strmatch(files(i,1),talker_group);  
    talker_number(i) = str2num(files(i,2:3)); 
end; 

%plot histogram
%{
%men
x = F0s(find(talker_group_code==1));   
figure(1);   
subplot(2,2,1);   
hist(x,20);  % use 20 “bins”   
set(gca,'XLim',[50 500]);  % set x-axis limits between 50 & 500 Hz   
%xlim([50 500]);
title('adult males')

%women
x = F0s(find(talker_group_code==2));
figure(1);
subplot(2,2,2);
hist(x,20);
set (gca, 'Xlim',[50 500]);
%xlim([50 500]);
title('adult females');

%boys
x = F0s(find(talker_group_code==3));
figure(1);
subplot(2,2,3);
hist(x,20);
set (gca, 'Xlim',[50 500]);
%xlim([50 500]);
title('boys');

%girls
x = F0s(find(talker_group_code==4));
figure(1);
subplot(2,2,4);
hist(x,20);
set (gca, 'Xlim',[50 500]);
%xlim([50 500]);
title('girls');
%}

%example calculate mean for all males
%{
x = F0s(find(talker_group_code==1)); 
mx = mean(x); 
disp('Mean F0 for males:') 
    disp(mx);
%}
    
%calculate the mean for each class set
%{
%ae
ae_tot = F0s(find(vowel_code==1));
m2 = mean(ae_tot(1:Ntrain));
%}

%finner middelverdien til f0 hver vokal, lagrer disse i en felles vektor
means_F0 = zeros(1,12);     %kan gjøre tilsvarende for F1, F2, F3
for i=1:12
    y = F0s(find(vowel_code==i));
    means_F0(i) = mean(y(1:Ntrain));
end

%middelverdivektor for F1
means_F1 = zeros(1,12);
for i=1:12
    y = F1s(find(vowel_code==i));
    means_F1(i) = mean(y(1:Ntrain));
end

%middelverdivektor for F2
means_F2 = zeros(1,12);
for i=1:12
    y = F2s(find(vowel_code==i));
    means_F2(i) = mean(y(1:Ntrain));
end


%prøver å lage en funksjon for å regne ut means
%{
function means = mean_tot(freq)
    means = zeros(1,12);     %kan gjøre tilsvarende for F1, F2, F3
    for i=1:12
        y = freq(find(vowel_code==i));
        means(i) = mean(y(1:Ntrain));
    end  
end
%}
