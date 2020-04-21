%husk diag() til diagonal matrise!!
clear;
N = 139;
Ntrain = 70;
Ntest = 69;
Nclass = 12;
Nfeatures = 3;
[files,dur,F0s,F1s,F2s,F3s,F4s,F120,F220,F320,F150,F250,F350,F180,F280,F380] = textread('vowdata_nohead.dat','%s%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f%4.1f'); 
nfiles = 1668;
Prob_w = 1/Nclass;
% character 1:     m=man, w=woman, b=boy, g=girl 
% characters 2-3:  talker number 
% characters 4-5:  vowel (ae=�had�, ah=�hod�, aw=�hawed�, eh=�head�,  
%    er=�heard�, ei=�haid�, ih=�hid�, iy=�heed�, oa=/o/ as in �boat�,  
%    oo=�hood�, uh=�hud�, uw=�who�d�) 

%databehandling
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
hist(x,20);  % use 20 �bins�   
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


%Oppgave 1 a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample mean for hver klasse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[mean_trainF1, mean_testF1] = find_mean(F1s,Nclass,N,Ntrain);
[mean_trainF2, mean_testF2] = find_mean(F2s,Nclass,N,Ntrain);
[mean_trainF3, mean_testF3] = find_mean(F3s,Nclass,N,Ntrain);


%12x3 matrise som inneholder gjennomsnittsverdien til f1, f2 og f3 for hver
%klasse
means_train = [mean_trainF1 mean_trainF2 mean_trainF3];
means_test = [mean_testF1 mean_testF2 mean_testF3];

%lager (12*70)x3 matrise for til trening av klassifisereren
Fs = make_string(F1s,F2s,F3s,N,Ntrain,Nclass);

% (12*3)x3 matrise best�ende av de 12 covarians matrisene
cov_matrices = zeros(Nclass*Nfeatures,Nfeatures);
cov_mat_test = zeros(Nclass*Nfeatures,Nfeatures);
for i = 1:Nclass
    cov_matrices((i-1)*3+1:(i*3),:) = find_cov(Fs, i, Ntrain);
    cov_mat_test((i-1)*3+1:(i*3),:) = find_cov(Fs, i, Ntest);
end

%cov_diag_train = cov_matrices;
%cov_diag_test = diag(cov_mat_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b: designe gaussian classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (70*12)x12 matrise som inneholder alle discriminant functions
g_all = zeros(Ntrain*Nclass,Nclass);
for i = 1:Ntrain*Nclass
    for c = 1:Nclass
        x = Fs(i,:);
        mu = means_train(c,:);
        cov_mat = cov_matrices((c-1)*Nfeatures+1:c*Nfeatures,:);
        g_all(i,c) = discriminant2(cov_mat, mu, x, Nfeatures, Prob_w);
    end
end



g_all_test = zeros(Ntest*Nclass, Nclass);
for i = 1:Ntest*Nclass
    for c = 1:Nclass
        x = Fs(i,:);
        mu = means_test(c,:);
        cov_mat = cov_matrices((c-1)*Nfeatures+1:c*Nfeatures,:);
        g_all_test(i,c) = discriminant2(cov_mat, mu, x, Nfeatures, Prob_w);
    end
end

%trener klassifisereren
true_val = fill_in_truevalues(Nclass,Ntrain);
trainset = test_classifier(Nclass, Ntrain,g_all,true_val);

%tester klassifisereren
true_val_test = fill_in_truevalues(Nclass,Ntest);
testset = test_classifier(Nclass, Ntest,g_all_test,true_val_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lage confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf_mat_train = compute_confusion(Nclass,Ntrain, true_val, trainset);
%{
%fill inn confusion matrix: 
conf_matrix= zeros(Nclass); % trenger en tabell som er 12x12, en med sann klasse, en med plassering. 
%length_test=length(testset_class); 
for t=1:Nclass*Ntrain
    x=true_val(t); 
    y=trainset(t); 
    conf_matrix(x,y)= conf_matrix(x,y) +1;
end
%}
disp(conf_mat_train);

conf_mat_test = compute_confusion(Nclass,Ntest, true_val_test, testset);
disp(conf_mat_test);
%finner error rate
error = compute_error(Nclass,Ntrain,conf_mat_train);
error_test = compute_error(Nclass,Ntest,conf_mat_test);
disp(error);
disp(error_test);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Funksjoner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%lager covariance matrix
function cov_matrix = find_cov(string, class_num, N)
    x_string = string((class_num*N-N)+1:class_num*N,:);
    cov_matrix = cov(x_string);
end

%funksjon 3.4 fra kompendiet
function g_i = discriminant2(cov_matrix, mu, x,Nfeatures, prior)
    frac = (sqrt(2*pi)^(Nfeatures)*det(cov_matrix))^(-1);
    expo = exp(-0.5*(x-mu)*cov_matrix^(-1)*(x-mu)');
    g_i = frac*expo*prior;
end

%Nclassx1 matrise best�ende av middelverdien til en feature for hver klasse
function [mean_train, mean_test] = find_mean(string,Nclass,Ntot,Ntrain)
    mean_train = zeros(Nclass,1);
    mean_test = zeros(Nclass,1);

    for i = 1:Nclass    
        y = string((i*Ntot-Ntot)+1:i*Ntot);
        mean_train(i,1) = mean(y(1:Ntrain));
        mean_test(i,1) = mean(y(Ntrain+1:Ntot));
    end
end

%matrise best�ende av inputs fra 3 ulike features for � gj�re behandlingen
%av data lettere
function matrix = make_string(string1,string2,string3,Ntot,Ntrain,Nclass)
    matrix_tot = [string1 string2 string3];
    matrix = zeros(Ntrain*Nclass,3);
    for i = 1:Nclass
        x_string = matrix_tot((i*Ntot-Ntot)+1:(i-1)*Ntot+Ntrain,:);
        matrix((i*Ntrain-Ntrain)+1:i*Ntrain,:) = x_string;
    end
end

%fyller inn sanne verdier til bruk til confusion matrix
function truevalues = fill_in_truevalues(C,N)
truevalues=zeros(C*N,1);
    for c = 1:C
        for i = 1:N
            truevalues(i+(c-1)*N,1) = c;
        end
    end
end

%lager confusion matrix
function conf_mat = compute_confusion(C,N, set_class, set_est)  
    conf_mat=zeros(C);
    for t=1:N*C
        x=set_class(t); 
        y=set_est(t); 
        conf_mat(x,y)=conf_mat(x,y)+1; 
    end
end

%finner error rate 
function error_rate = compute_error(C,N,conf_mat)
sum_error=0;
for i = 1:C
    for j = 1:C
        if j ~= i           
            sum_error = sum_error +conf_mat(j,i);
        end
    end
end
error_rate = sum_error/(C*N);
end

%trener klassifisereren, fyller inn resultatene fra discriminant functions
function testset = test_classifier(C,N,g_all,true_set)
testset = zeros(N*C,1);
for x = 1:N*C
    for c = 1:C
        if g_all(x,c) == max(g_all(x,:))
            testset(x,1)= c;
        end
    end
end
end

%{
function g_all = make_gvec(C,N,features,means,cov_mats,discriminant2)
g_all = zeros(N*C,C);
for i = 1:N*C
    for c = 1:C
        x = features(i,:);
        mu = means(c,:);
        cov_mat = cov_mats((c-1)*N+1:c*N,:);
        g_all(i,c) = discriminant2(cov_mat, mu, x, Nfeatures, Prob_w);
    end
end
end
%}