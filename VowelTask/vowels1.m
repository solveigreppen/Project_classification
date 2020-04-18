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
% characters 4-5:  vowel (ae=”had”, ah=”hod”, aw=”hawed”, eh=”head”,  
%    er=”heard”, ei=”haid”, ih=”hid”, iy=”heed”, oa=/o/ as in “boat”,  
%    oo=”hood”, uh=”hud”, uw=”who’d”) 


vowel = str2mat('ae','ah','aw','eh','er','ei','ih','iy','oa','oo','uh','uw'); 
talker_group = char('m','w','b','g');  

vowel_code = zeros(1,nfiles);
talker_number = zeros(1,nfiles);

files=char(files); % convert cell array to character matrix [nfiles,nchar]=size(filenames); 
for i=1:nfiles  
    vowel_code(i) = strmatch(files(i,4:5),vowel);  %legger inn en nummerkode for hver av vowelene i en fil
    talker_group_code(i) = strmatch(files(i,1),talker_group);  %legger inn en nummerkode avhengig av hvem som prater
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


%Oppgave 1 a)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sample mean for hver klasse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

<<<<<<< Updated upstream
[mean_trainF150, mean_test150] = find_mean(F150,Nclass,N,Ntrain);
[mean_trainF250, mean_test250] = find_mean(F250,Nclass,N,Ntrain);
[mean_trainF350, mean_test350] = find_mean(F350,Nclass,N,Ntrain);

%{
%finner middelverdien til f1 for hver vokalklasse, lagrer disse i en felles vektor
=======
%finner middelverdien til f1 for hver vokalklasse, lagrer disse i en felles
%vektor, f1 er feature 1 osv til feature 3
>>>>>>> Stashed changes
means_F1_train = zeros(Nclass,1);     %kan gjøre tilsvarende for F1, F2, F3
means_F1_test = zeros(Nclass,1);
for i=1:12
    y = F1s(find(vowel_code==i));
    means_F1_train(i,1) = mean(y(1:Ntrain));
    means_F1_test(i,1) = mean(y(Ntrain+1:N));
end

<<<<<<< Updated upstream
%middelverdivektor for F2
means_F2_train = zeros(Nclass,1);     
=======
%middelverdivektor for F2 - kan kanskje lage funksjon
means_F2_train = zeros(Nclass,1);     %kan gjøre tilsvarende for F1, F2, F3
>>>>>>> Stashed changes
means_F2_test = zeros(Nclass,1);
for i=1:12
    y = F2s(find(vowel_code==i));
    means_F2_train(i,1) = mean(y(1:Ntrain));
    means_F2_test(i,1) = mean(y(Ntrain+1:N));
end

%middelverdivektor for F3
means_F3_train = zeros(Nclass,1);     
means_F3_test = zeros(Nclass,1);
for i=1:12
    y = F3s(find(vowel_code==i));
    means_F3_train(i,1) = mean(y(1:Ntrain));
    means_F3_test(i,1) = mean(y(Ntrain+1:N));
end
%}

%12x3 matrise som inneholder gjennomsnittsverdien til f1, f2 og f3 for hver
%klasse
means_train = [mean_trainF150 mean_trainF250 mean_trainF350];
means_test = [mean_test150 mean_test250 mean_test350];

%lager (12*70)x3 matrise for til trening av klassifisereren
f_50 = make_string(F150,F250,F350,N,Ntrain,Nclass);

%{
Fs = [F1s F2s F3s];

F_50 = [F150 F250 F350];
F_50_train = zeros(Ntrain*Nclass,3);
for i = 1:Nclass
    x_string = F_50((i*N-N)+1:(i-1)*N+Ntrain,:);
    F_50_train((i*Ntrain-Ntrain)+1:i*Ntrain,:) = x_string;
end

%lager (12*70)x3 matrise for trainset
strings = zeros(70*Nclass,3);
for i = 1:Nclass
    x_string = Fs((i*N-N)+1:(i-1)*N+Ntrain,:);
    strings((i*Ntrain-Ntrain)+1:i*Ntrain,:) = x_string;
end
%}
% (12*3)x3 matrise bestående av de 12 covarians matrisene
cov_matrices = zeros(Nclass*Nfeatures,Nfeatures);
for i = 1:Nclass
    cov_matrices((i-1)*3+1:(i*3),:) = find_cov(f_50, i, Ntrain);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1b: designe gaussian classifier
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g_all = zeros(3*12,3);
for i = 1:Nclass
    g_all((i*3-3)+1:i*3,:) = discriminant2(cov_matrices((i*3-3)+1:i*3,:),means_train(i,:)',f_50(1,:),Prob_w);
end
%testing the functions
%{
g_840_1 = discriminant(Nfeatures, cov_matrices(1:3,:), means_train(1,:)',strings(840,:));
g_840_12 = discriminant(Nfeatures, cov_matrices(34:36,:), means_train(12,:)',strings(840,:));
g_1 = discriminant(Nfeatures,cov_matrices(1:3,1:3),means_train(1,:)',strings(1,:));
g_1_1 = discriminant2(cov_matrices(1:3,1:3),means_train(1,:)',strings(1,:));
g_1_12 = discriminant2(cov_matrices(34:36,1:3),means_train(12,:)',strings(1,:));
g_1_8 = discriminant2(cov_matrices(22:24,1:3),means_train(8,:)',strings(1,:));
%}
%{
gs = zeros (3*12,3);
for i = 1:Nclass
    gs((i*3-3)+1:i*3,:) = discriminant2(cov_matrices((i*3-3)+1:i*3,:),means_train(i,:)',strings(1,:))*Prob_w;
end
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lage confusion matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%kode fra iris
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

%kode fra iris
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

%{
mean_ae = find_mean(1,Ntrain,F1s);
function mean_i = find_mean(vowel_num, N, feature)
    y = feature(find(vowel_code==vowel_num));
    mean_i = mean(y(1:N));
end
%}
%{
function x_string = find_x(Ntot, Ntrain, class_num,vec)
    x_string = vec((class_num*Ntot-Ntrain)+1:(class_num-1)*Ntot+Ntrain,:);
end
%}
function cov_matrix = find_cov(string, class_num, N)
    x_string = string((class_num*N-N)+1:class_num*N,:);
    cov_matrix = cov(x_string);
end

function g_i = discriminant(dim,cov_matrix, mu,x)
    g_i = -(dim/2)*log(2*pi)-0.5*log(abs(cov_matrix))-0.5*(x-mu)'*cov_matrix^(-1)*(x-mu);
end

function g_i = discriminant2(cov_matrix, mu, x, prior)
    g_i = normpdf(x,mu,cov_matrix)*prior;
end


function [mean_train, mean_test] = find_mean(string,Nclass,Ntot,Ntrain)
    mean_train = zeros(Nclass,1);
    mean_test = zeros(Nclass,1);

    for i = 1:Nclass    
        y = string((i*Ntot-Ntot)+1:i*Ntot);
        mean_train(i,1) = mean(y(1:Ntrain));
        mean_test(i,1) = mean(y(Ntrain+1:Ntot));
    end
end

function matrix = make_string(string1,string2,string3,Ntot,Ntrain,Nclass)
    matrix_tot = [string1 string2 string3];
    matrix = zeros(Ntrain*Nclass,3);
    for i = 1:Nclass
        x_string = matrix_tot((i*Ntot-Ntot)+1:(i-1)*Ntot+Ntrain,:);
        matrix((i*Ntrain-Ntrain)+1:i*Ntrain,:) = x_string;
    end
end

function conf_mat = compute_confusion(N, C, set_class, set_est) %trengs til oppg 1b
conf_mat=zeros(C);
    for t=1:N*C
        x=set_class(t); 
        y=set_est(t); 
        conf_mat(x,y)=conf_mat(x,y)+1; 
    end
end

function error = error_rate(C, N, conf_mat) %trengs til opp 1b)
sum_error=0;
for i = 1:C
    for j = 1:C
        if j ~= i           
            sum_error = sum_error +conf_mat(j,i);
        end
    end
end
error = sum_error/(C*N);
end

