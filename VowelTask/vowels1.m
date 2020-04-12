%husk diag() til diagonal matrise!!
clear;
D = 139;
Ntrain = 70;
Ntest = 69;
Nclass = 12;
Nfeatures = 3;
[files,dur,F0s,F1s,F2s,F3s,F4s,F120,F220,F320,F150,F250,F350,F180,F280,F380] = textread('vowdata_nohead.dat','%s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f'); 
nfiles = 1668;
Prob_w = 1/Nclass;
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

%(12*139)x12 matrise.
vowels = zeros(D*Nclass,Nclass);
for i = 1:Nclass
    vowels(1:D,i) = F1s(find(vowel_code==i));
    vowels(D+1:2*D,i) = F2s(find(vowel_code==i));
    vowels(2*D+1:3*D,i) = F3s(find(vowel_code==i));
end

%train_vowels = zeros()
%test_vowels = vowels(Ntrain+1:D,:);


%Oppgave 1 a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample mean for hver klasse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%finner middelverdien til f1 for hver vokalklasse, lagrer disse i en felles vektor
means_F1 = zeros(Nclass,1);     %kan gjøre tilsvarende for F1, F2, F3
for i=1:12
    y = F1s(find(vowel_code==i));
    means_F1(i,1) = mean(y(1:Ntrain));
end

%middelverdivektor for F2
means_F2 = zeros(Nclass,1);
for i=1:12
    y = F2s(find(vowel_code==i));
    means_F2(i,1) = mean(y(1:Ntrain));
end

%middelverdivektor for F3
means_F3 = zeros(Nclass,1);
for i=1:12
    y = F3s(find(vowel_code==i));
    means_F3(i,1) = mean(y(1:Ntrain));
end

%12x3 matrise som inneholder gjennomsnittsverdien til f1, f2 og f3 for hver
%klasse
means = [means_F1 means_F2 means_F3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lager confusion matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
cov_matrices = zeros(Nclass*Nfeatures,Nfeatures); %matrise som inneholder alle 12 cov.matrisene
for c = 1:Nclass %1-12
    for i = 1:Nfeatures %1-3
        u_i = means(c,i);
        for k = 1:Ntrain
            x_k = 
        end
        
    end
end
%}
%{
cov_matrices = zeros(Nclass*Nfeatures,Nfeatures); %matrise som inneholder alle 12 cov.matrisene
for c = 1:Nfeatures
    cov_matrix = zeros(Nfeatures);
    for k = 1:Ntrain
        for i = 1:Nclass
            u_i = means(i,c);
            x_k = F1s(k);
            cov_matrix = cov_matrix + (x_k-u_i)*((x_k-u_i).');
        end
    end
    cov_matrix = cov_matrix./Ntrain;
    cov_matrices((i-1)*3+1:(i*3),:) = cov_matrix;
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b: designe gaussian classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
means_test = zeros(1,Nclass);
for i=1:12
    y = F1s(find(vowel_code==i));
    means_test(i) = mean(y(Ntrain+1:D));
end

prob_x_w = zeros(Nfeatures*Nclass,Nfeatures);
for i = 1:Nclass
    u_i = means_F1(1,i);
    x = F1s(1);
    sigma = cov_matrices((i-1)*Nfeatures+1:i*Nfeatures,:);
    prob = normpdf(x,u_i,sigma);
    prob_x_w((i-1)*Nfeatures+1:i*Nfeatures,:) = prob;
end
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lage confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
g_all_test = zeros(1,Nclass);
trainset_class = zeros(1,Nclass);
testset_est = zeros(1,Nclass);

for x = 1:C*Ntest
    for c = 1:C
        if g_all_test(x,c) == max(g_all_test(x,:))
            testset_est(x)= c;
        end
    end
end
%}

%{
%cov.matrix for test set 
cov_matrices_test = zeros(Nfeatures,Nclass*Nfeatures); %matrise som inneholder alle 12 cov.matrisene
 %matrise som inneholder alle 12 cov.matrisene
for i = 1:Nclass
    cov_matrix = zeros(Nfeatures);
    for k = Ntrain+1:D
        u_i = means_test(i);
        x_k = F1s(k);
        cov_matrix = cov_matrix + (x_k-u_i)*((x_k-u_i).');
    end
    cov_matrix = cov_matrix./Ntest;
    cov_matrices_test(:,(i-1)*3+1:(i*3)) = cov_matrix;
end
C1 = cov_matrices(1:3,1:3);
C2 = cov_matrices_test(1:3,1:3);
%C = confusionmat(C1,C2);
%}
%{
conf_matrix= zeros(Nclass); % trenger en tabell som er 12x12, en med sann klasse, en med plassering. 
%fill inn confusion matrix: 
%length_test=length(testset_class); 
for t=1:length_test
    x=testset_class(t); 
    y=testset_est(t); 
    conf_matrix(x,y)= conf_matrix(x,y) +1;
end
disp(conf_matrix);
%}
